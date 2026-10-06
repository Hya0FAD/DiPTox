# diptox/core.py
import os
import sys

import pandas as pd
from typing import Optional, List, Union, Tuple, Callable, Dict, Any
from functools import partial, wraps
from copy import deepcopy
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import Descriptors
from datetime import datetime
import requests
import multiprocessing as mp
import platform
import re
from .chem_processor import ChemistryProcessor, AmbiguousFragmentError
from .element_policy import is_metal
from .web_request import WebService
from .data_io import DataHandler
from .data_deduplicator import DataDeduplicator
from .substructure_search import SubstructureSearcher
from .unit_processor import UnitProcessor, TARGET_SCALES_ATTR, get_target_scale
from .column_policy import INITIAL_COLUMNS, PREPROCESS_COLUMNS, source_column_renames
from .logger import log_manager
logger = log_manager.get_logger(__name__)
from diptox import user_reg


def _step_failure_message(description: str, mol: Chem.Mol) -> str:
    """Return a record-level explanation for a failed standardization step."""
    if description == "Dummy atom check":
        dummy_count = sum(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms())
        return f"Dummy atom check failed: structure contains {dummy_count} wildcard atom(s)"
    if description == "Mixture removal":
        fragment_count = len(Chem.GetMolFrags(mol))
        return f"Mixture removal failed: {fragment_count} unresolved components"
    if description == "Metal check":
        metals = sorted({
            atom.GetSymbol() for atom in mol.GetAtoms()
            if is_metal(atom.GetAtomicNum())
        })
        return f"Metal check failed: metal-containing structure ({', '.join(metals)})"
    if description == "Inorganic removal":
        if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
            return "Inorganic removal failed: structure contains no carbon"
        return "Inorganic removal failed: carbon-containing inorganic exclusion rule matched"
    if description == "Radical check":
        radical_atoms = sorted({
            atom.GetSymbol() for atom in mol.GetAtoms()
            if atom.GetNumRadicalElectrons() > 0
        })
        return f"Radical check failed: unsupported radical center ({', '.join(radical_atoms)})"
    if description == "Non-neutral rejection":
        return f"Non-neutral rejection failed: total formal charge is {Chem.GetFormalCharge(mol):+d}"
    return f"{description} failed"


def _standardization_status(is_valid: bool, processing_log: str) -> str:
    if is_valid:
        return "Retained"
    status_markers = [
        ("Ambiguous parent", "Ambiguous parent"),
        ("Empty SMILES", "Empty input"),
        ("Invalid SMILES", "Invalid structure"),
        ("Dummy atom check", "Dummy atom excluded"),
        ("Mixture removal", "Mixture excluded"),
        ("Metal check", "Metal excluded"),
        ("Inorganic removal", "Inorganic excluded"),
        ("Radical check", "Radical excluded"),
        ("Non-neutral rejection", "Charge excluded"),
        ("Atom validation", "Element excluded"),
    ]
    for marker, status in status_markers:
        if marker in processing_log:
            return status
    return "Processing error"


def _molecule_metadata(smiles: Any, sanitize: bool = True) -> Dict[str, Any]:
    """Calculate stable record-level metadata used for parent-mapping audits."""
    if not isinstance(smiles, str) or not smiles.strip():
        return {}
    try:
        parser = Chem.SmilesParserParams()
        parser.sanitize = sanitize
        parser.removeHs = False
        mol = Chem.MolFromSmiles(smiles, parser)
        if mol is None:
            return {}
        if not sanitize:
            mol.UpdatePropertyCache(strict=False)
        metals = sorted({
            atom.GetSymbol() for atom in mol.GetAtoms()
            if is_metal(atom.GetAtomicNum())
        })
        return {
            'smiles': Chem.MolToSmiles(mol, canonical=True),
            'identity': ChemistryProcessor.standardize_smiles(mol),
            'fragments': len(Chem.GetMolFrags(mol)),
            'charge': Chem.GetFormalCharge(mol),
            'mw': Descriptors.MolWt(mol),
            'metals': metals,
        }
    except Exception:
        return {}


def check_data_loaded(func):
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        if platform.system() == "Windows" and mp.current_process().name != 'MainProcess':
            return func(self, *args, **kwargs)

        if self.df is None:
            raise ValueError("No data loaded. Please load data first.")
        return func(self, *args, **kwargs)
    return wrapper


def _run_on_main_process_only(func):
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        if platform.system() == "Windows" and mp.current_process().name != 'MainProcess':
            return None
        return func(self, *args, **kwargs)
    return wrapper


def _worker_preprocess(args: Tuple[str, ChemistryProcessor, Dict[str, Any]]) -> Tuple[bool, str, Optional[str]]:
    """
    Independent worker function for parallel processing.
    Must be defined at module level to be picklable by multiprocessing.
    """
    try:
        from rdkit import RDLogger
        RDLogger.DisableLog('rdApp.*')
    except (ImportError, AttributeError):
        pass
    smiles, processor, config = args
    comments = []

    # Validate one cell without applying scalar null checks to lists/arrays.
    if not isinstance(smiles, str):
        if pd.api.types.is_scalar(smiles) and pd.isna(smiles):
            return False, "Empty SMILES", None
        return False, "Invalid SMILES type: expected a string", None
    if not smiles.strip():
        return False, "Empty SMILES", None

    try:
        # 1. Initialization
        mol = processor.smiles_to_mol(smiles, config['sanitize'], remove_hs=False)
        if mol is None:
            return False, "Invalid SMILES", None

        # 2. Build pipeline dynamically based on config
        # Note: We reconstruct the steps here because passing partial objects
        # across processes can sometimes be problematic or inefficient.
        steps = []
        step_descriptions = []

        steps.append(processor.reject_dummy_atoms)
        step_descriptions.append("Dummy atom check")

        if config['remove_stereo']:
            steps.append(processor.remove_stereochemistry)
            step_descriptions.append("Stereo removal")

        if config['remove_isotopes']:
            steps.append(processor.remove_isotopes)
            step_descriptions.append("Isotope removal")

        if config['remove_hs']:
            steps.append(processor.remove_hydrogens)
            step_descriptions.append("Hydrogen removal")

        if config['remove_salts']:
            steps.append(processor.remove_salts)
            step_descriptions.append("Salt removal")

        if config['remove_solvents']:
            steps.append(processor.remove_solvents)
            step_descriptions.append("Solvent removal")

        worker_mixture_mode = config.get('mixture_mode')
        if worker_mixture_mode is None:
            worker_mixture_mode = (
                "keep" if not config.get('remove_mixtures', True)
                else ("largest" if config.get('keep_largest_fragment', False) else "reject")
            )
        mixture_processing_enabled = (
            config.get('remove_mixtures', True) or worker_mixture_mode != "keep"
        )
        if mixture_processing_enabled:
            mixture_processor = partial(
                processor.remove_mixtures,
                hac_threshold=config['hac_threshold'],
                mode=worker_mixture_mode,
            )
            steps.append(mixture_processor)
            step_descriptions.append(
                "Mixture preservation" if worker_mixture_mode == "keep" else "Mixture removal"
            )
        else:
            steps.append(lambda mol: mol)
            step_descriptions.append("Mixture preservation")

        if config['reject_metal_species']:
            steps.append(processor.reject_metals)
            step_descriptions.append("Metal check")

        if config['remove_inorganic']:
            steps.append(processor.remove_inorganic)
            step_descriptions.append("Inorganic removal")

        if config['reject_radical_species']:
            steps.append(processor.reject_radicals)
            step_descriptions.append("Radical check")

        if config['neutralize']:
            charge_processor = partial(
                processor.neutralize_charges,
                reject_non_neutral=False
            )
            steps.append(charge_processor)
            step_descriptions.append("Charge neutralization")

        if config['reject_non_neutral']:
            steps.append(processor.reject_non_neutral)
            step_descriptions.append("Non-neutral rejection")

        if config['check_valid_atoms']:
            steps.append(processor.effective_atom)
            step_descriptions.append("Atom validation")

        if config.get('add_hs', False):
            steps.append(processor.add_hydrogens)
            step_descriptions.append("Hydrogen addition")

        if config['remove_salts'] or config['remove_solvents'] or mixture_processing_enabled:
            steps.append(processor.collapse_identical_components)
            step_descriptions.append("Identical component collapse")

        # 3. Execute Pipeline
        for step, desc in zip(steps, step_descriptions):
            if mol is None:
                break
            try:
                before = processor.standardize_smiles(mol)
                processed = step(mol)
                if processed is None:
                    comments.append(_step_failure_message(desc, mol))
                    mol = None
                else:
                    mol = processed
                    after = processor.standardize_smiles(mol)
                    if desc == "Mixture preservation" and len(Chem.GetMolFrags(mol)) > 1:
                        comments.append(
                            f"Mixture preservation: retained {len(Chem.GetMolFrags(mol))} components"
                        )
                    elif before != after:
                        comments.append(f"{desc}: {before} -> {after}")
            except AmbiguousFragmentError as e:
                comments.append(f"Ambiguous parent: {e}")
                mol = None
            except Exception as e:
                comments.append(f"{desc} error: {str(e)}")
                mol = None

        # 4. Finalize
        if mol is not None:
            try:
                std_smiles = processor.standardize_smiles(mol)
                mapped_atoms = sum(atom.GetAtomMapNum() != 0 for atom in mol.GetAtoms())
                if mapped_atoms:
                    comments.append(
                        f"Atom map removal: omitted {mapped_atoms} mapping label(s) from Canonical SMILES"
                    )
                return True, "; ".join(comments) if comments else "Success", std_smiles
            except Exception as e:
                return False, f"Standardization failed: {str(e)}", None
        else:
            return False, "; ".join(comments), None

    except Exception as e:
        return False, f"Worker Error: {str(e)}", None


class DiptoxPipeline:
    """Main processing class that coordinates various modules."""

    def __init__(self, interactive: bool = True):
        """Initialize the pipeline; disable interactive prompts for automation."""
        self.interactive = interactive
        if not interactive or (platform.system() == "Windows" and mp.current_process().name != 'MainProcess'):
            pass
        else:
            self._check_initial_registration()
        self.chem_processor = ChemistryProcessor()
        self.data_handler = DataHandler()

        self.df: Optional[pd.DataFrame] = None
        self.source_columns: List[str] = []
        self.excluded_df: pd.DataFrame = pd.DataFrame()
        self.smiles_col: str = "Smiles"
        self.cas_col: Optional[str] = None
        self.name_col: Optional[str] = None
        self.target_col: Optional[str] = None
        self.unit_col: Optional[str] = None
        self.inchikey_col = None
        self.id_col: Optional[str] = None

        self.deduplicator = None
        self.web_service = None
        self._preprocess_key = 0
        self._units_standardized = False
        self._dedup_unit_settings = None
        self.web_source = None
        self._audit_log = []
        self._history = []
        self._max_history = 5
        self._current_dedup_config = None
        self._structure_derived_columns = set()
        # Stable row identities live only in the index, never in exported columns.
        self._source_index_labels = {}
        self._row_lineage = {}
        self._input_column_aliases = {}

    @staticmethod
    def _check_initial_registration():
        """
        Check if user info is registered.
        Only triggers in interactive terminal mode AND not in GUI mode.
        """
        if os.environ.get("DIPTOX_GUI_MODE") == "true":
            return

        if user_reg.is_registered_or_skipped():
            return

        try:
            print("=" * 60)
            print("👋 Welcome to DiPTox!")
            print("To help us improve, we'd love to know who is using this tool.")
            print("This is OPTIONAL. You can skip by pressing Enter directly.")
            print("=" * 60)

            name = input("Enter your Name (or press Enter to skip): ").strip()
            if not name:
                print("Skipping registration. You won't be asked again.")
                user_reg.save_status("skipped")
                print("=" * 60)
                return

            affiliation = input("Enter your Affiliation/Unit: ").strip()
            email = input("Enter your Email (optional): ").strip()

            print("Sending information...", end=" ")
            success, msg = user_reg.submit_info(name, affiliation, email)
            if success:
                print("Done! Thank you.")
            else:
                print(f"(Note: {msg})")
            print("=" * 60 + "\n")

        except (EOFError, OSError):
            pass

    def _record_step(self, step_name: str, df_before: Optional[pd.DataFrame], df_after: pd.DataFrame,
                     details: str = ""):
        """Records a processing step into the audit log."""
        rows_before = len(df_before) if df_before is not None else 0
        rows_after = len(df_after) if df_after is not None else 0

        # Calculate delta. For initial load, before is 0, so delta is +rows.
        # For filtering steps, delta is negative (rows removed).
        if df_before is None:
            delta = f"+{rows_after}"
        else:
            diff = rows_after - rows_before
            delta = str(diff) if diff <= 0 else f"+{diff}"

        entry = {
            "Step": step_name,
            "Timestamp": datetime.now().strftime("%H:%M:%S"),
            "Rows Before": rows_before,
            "Rows After": rows_after,
            "Delta": delta,
            "Details": details
        }
        self._audit_log.append(entry)

    def _save_checkpoint(self):
        """A snapshot of the current data and the status of key columns"""
        if self.df is not None:
            snapshot = {
                'df': self.df.copy(),
                'excluded_df': self.excluded_df.copy(),
                'smiles_col': self.smiles_col,
                'cas_col': self.cas_col,
                'name_col': self.name_col,
                'target_col': self.target_col,
                'unit_col': self.unit_col,
                'id_col': self.id_col,
                'inchikey_col': self.inchikey_col,
                '_preprocess_key': self._preprocess_key,
                '_units_standardized': self._units_standardized,
                '_structure_derived_columns': set(self._structure_derived_columns),
                '_source_index_labels': dict(self._source_index_labels),
                '_row_lineage': dict(self._row_lineage),
                'source_columns': list(self.source_columns),
                '_input_column_aliases': dict(self._input_column_aliases),
                '_dedup_unit_settings': deepcopy(self._dedup_unit_settings),
                'deduplicator': deepcopy(self.deduplicator),
            }
            self._history.append(snapshot)
            if len(self._history) > self._max_history:
                self._history.pop(0)

    @_run_on_main_process_only
    def undo(self) -> bool:
        if not hasattr(self, '_history') or not self._history:
            return False

        df_before_undo = self.df.copy() if self.df is not None else None
        snapshot = self._history.pop()
        self.df = snapshot['df']
        self.excluded_df = snapshot.get('excluded_df', pd.DataFrame()).copy()
        self.smiles_col = snapshot['smiles_col']
        self.cas_col = snapshot['cas_col']
        self.name_col = snapshot['name_col']
        self.target_col = snapshot['target_col']
        self.unit_col = snapshot['unit_col']
        self.id_col = snapshot['id_col']
        self.inchikey_col = snapshot.get('inchikey_col')
        self._preprocess_key = snapshot['_preprocess_key']
        self._units_standardized = snapshot['_units_standardized']
        self._structure_derived_columns = snapshot.get('_structure_derived_columns', set())
        self._source_index_labels = snapshot.get('_source_index_labels', {}).copy()
        self._row_lineage = snapshot.get('_row_lineage', {}).copy()
        self.source_columns = snapshot.get('source_columns', self.source_columns).copy()
        self._input_column_aliases = snapshot.get('_input_column_aliases', {}).copy()
        self._dedup_unit_settings = deepcopy(snapshot.get('_dedup_unit_settings'))
        self.deduplicator = deepcopy(snapshot.get('deduplicator', self.deduplicator))

        self._record_step("Undo", df_before_undo, self.df, "Restored to previous state")
        return True

    @_run_on_main_process_only
    def get_history(self) -> pd.DataFrame:
        """Returns the processing history as a DataFrame."""
        if not self._audit_log:
            return pd.DataFrame(columns=["Step", "Timestamp", "Rows Before", "Rows After", "Delta", "Details"])
        return pd.DataFrame(self._audit_log)

    @_run_on_main_process_only
    def load_data(self,
                  input_data: Union[str, List[str], pd.DataFrame],
                  smiles_col: str = None,
                  cas_col: Optional[str] = None,
                  name_col: Optional[str] = None,
                  target_col: Optional[str] = None,
                  unit_col: Optional[str] = None,
                  inchikey_col: Optional[str] = None,
                  id_col: Optional[str] = None,
                  **kwargs) -> None:
        """
        Load data and initialize columns for processing.
        :param input_data: Path to input data (.xlsx/.xls/.csv/.txt/.sdf/.smi), or a list, or a DataFrame.
        :param smiles_col: The column name containing SMILES strings (optional).
        :param cas_col: The column name for CAS Numbers (optional).
        :param name_col: The column name for names (optional).
        :param target_col: The column name for target values (optional).
        :param unit_col: The column name for units of the target values (optional).
        :param inchikey_col: The column name for InChIKeys (optional).
        :param id_col: The column name for SMI file's SMILES ID (optional)
        :param sep: CSV file delimiter.
        """
        if smiles_col is None and isinstance(input_data, str) and input_data.lower().endswith(('.sdf', '.smi')):
            smiles_col = 'smiles'

        user_specified_smiles = smiles_col
        df = self.data_handler.load_data(input_data=input_data, smiles_col=smiles_col, cas_col=cas_col,
                                         name_col=name_col,
                                         target_col=target_col, unit_col=unit_col, inchikey_col=inchikey_col,
                                         id_col=id_col, interactive=self.interactive, **kwargs)

        self.source_columns = list(df.columns)
        self._source_index_labels = dict(enumerate(df.index.tolist()))
        df.index = pd.RangeIndex(len(df))
        self._row_lineage = {index: (index,) for index in df.index}
        initial_renames = source_column_renames(df.columns, self.source_columns, INITIAL_COLUMNS)
        df.rename(columns=initial_renames, inplace=True)
        self.source_columns = [initial_renames.get(column, column) for column in self.source_columns]
        self._input_column_aliases = dict(initial_renames)
        user_specified_smiles = initial_renames.get(user_specified_smiles, user_specified_smiles)
        cas_col = initial_renames.get(cas_col, cas_col)
        name_col = initial_renames.get(name_col, name_col)
        target_col = initial_renames.get(target_col, target_col)
        unit_col = initial_renames.get(unit_col, unit_col)
        inchikey_col = initial_renames.get(inchikey_col, inchikey_col)
        id_col = initial_renames.get(id_col, id_col)
        if 'Canonical SMILES' not in df.columns:
            df['Canonical SMILES'] = pd.Series(pd.NA, index=df.index, dtype="string")
        if 'Processing Log' not in df.columns:
            df['Processing Log'] = pd.Series(pd.NA, index=df.index, dtype="string")
        if 'Standardization Status' not in df.columns:
            df['Standardization Status'] = pd.Series(pd.NA, index=df.index, dtype="string")
        if 'Is Valid' not in df.columns:
            df['Is Valid'] = pd.Series(False, index=df.index, dtype="boolean")

        try:
            df['Canonical SMILES'] = df['Canonical SMILES'].astype("string")
        except Exception:
            df['Canonical SMILES'] = pd.Series(pd.NA, index=df.index, dtype="string")
        try:
            df['Processing Log'] = df['Processing Log'].astype("string")
        except Exception:
            df['Processing Log'] = pd.Series(pd.NA, index=df.index, dtype="string")
        try:
            df['Standardization Status'] = df['Standardization Status'].astype("string")
        except Exception:
            df['Standardization Status'] = pd.Series(pd.NA, index=df.index, dtype="string")
        try:
            df['Is Valid'] = df['Is Valid'].astype("boolean")
        except Exception:
            df['Is Valid'] = pd.Series(False, index=df.index, dtype="boolean")

        self.df = df
        self.excluded_df = pd.DataFrame()
        if user_specified_smiles:
            self.smiles_col = user_specified_smiles
        else:
            if isinstance(input_data, str) and input_data.lower().endswith('.sdf'):
                self.smiles_col = 'smiles'
                logger.info("SDF file detected: 'smiles' column automatically assigned.")
            else:
                self.smiles_col = None
        self.cas_col = cas_col
        self.name_col = name_col
        self.target_col = target_col
        self.unit_col = unit_col
        self.inchikey_col = inchikey_col
        self.id_col = id_col
        self._preprocess_key = 0
        self._units_standardized = False
        self._dedup_unit_settings = None
        self.deduplicator = None
        self._current_dedup_config = None
        self._structure_derived_columns = set()

        if hasattr(self, '_history'):
            self._history.clear()
        else:
            self._history = []
        self._audit_log = []

        source_name = input_data if isinstance(input_data, str) else "Memory/List"
        self._record_step("Data Loading", None, self.df, f"Source: {source_name}")
        if 'Structure Reliability' in self.df:
            unreliable = self.df['Structure Reliability'].eq('Unreliable').fillna(False)
            if unreliable.any():
                excluded = self.df.loc[unreliable].copy()
                excluded['Is Valid'] = False
                excluded['Standardization Status'] = 'Unreliable structure'
                reasons = excluded.get('Import Error', pd.Series(pd.NA, index=excluded.index))
                excluded['Processing Log'] = reasons.fillna('SDF structure and declared SMILES disagree')
                excluded['Exclusion Reason'] = excluded['Processing Log']
                excluded['Excluded At'] = 'Input validation'
                excluded['Source DataFrame Index'] = self._source_labels(excluded.index)
                self._merge_excluded_records(excluded)
                df_before = self.df
                self.df = self.df.loc[~unreliable].copy()
                self._record_step(
                    'Input validation', df_before, self.df,
                    f'Excluded {len(excluded)} unreliable SDF structure record(s)',
                )

    def _source_labels(self, indices) -> List[Any]:
        """Translate private row identities to the original input index labels."""
        return [self._source_index_labels.get(index, index) for index in indices]

    def _protect_source_columns(self, output_columns):
        """Keep imported values and update mappings before writing derived fields."""
        renamed = source_column_renames(self.df.columns, self.source_columns, output_columns)
        if not renamed:
            return {}
        self.df.rename(columns=renamed, inplace=True)
        if not self.excluded_df.empty:
            self.excluded_df.rename(columns=renamed, inplace=True)
        self.source_columns = [renamed.get(column, column) for column in self.source_columns]
        self._input_column_aliases = {
            key: renamed.get(value, value) for key, value in self._input_column_aliases.items()
        }
        self._input_column_aliases.update(renamed)
        for role in ('smiles_col', 'target_col', 'unit_col', 'cas_col', 'name_col', 'id_col', 'inchikey_col'):
            setattr(self, role, renamed.get(getattr(self, role), getattr(self, role)))
        if self.deduplicator:
            self.deduplicator.condition_cols = [
                renamed.get(column, column) for column in self.deduplicator.condition_cols
            ]
        return renamed

    def _invalidate_standardized_mw_targets(self):
        """Clear derived values whose MW basis no longer matches the current structure."""
        dependencies = self.df.attrs.get('_diptox_standardized_mw_dependencies', {})
        stale_targets = set(self.df.attrs.get('_diptox_stale_targets', []))
        for column, weights in dependencies.items():
            if column not in self.df:
                continue
            stale = [index for index in self.df.index if index in weights and (
                pd.isna(self.df.at[index, 'Standardized Molecular Weight'])
                or abs(float(self.df.at[index, 'Standardized Molecular Weight']) - weights[index]) > 1e-6
            )]
            if stale:
                self.df.loc[stale, column] = pd.NA
                stale_targets.add(column)
                if column == self.target_col:
                    self.df.loc[stale, 'Unit Conversion Status'] = 'Stale standardized molecular weight'
                    self._units_standardized = False
        self.df.attrs['_diptox_stale_targets'] = sorted(stale_targets)

    def _restore_target_scale_from_units(self):
        """Recover an exported logarithmic target's scale without extra user mapping."""
        if not self.target_col or not self.unit_col or self.unit_col not in self.df:
            return
        if self.target_col in self.df.attrs.get(TARGET_SCALES_ATTR, {}):
            return
        units = self.df[self.unit_col].dropna().astype(str)
        scales = set()
        for unit in units:
            match = re.fullmatch(r'(-?log10)\(.*\)', unit)
            scales.add(match.group(1) if match else 'linear')
        if len(scales) > 1 and scales != {'linear'}:
            raise ValueError("The target mixes linear and logarithmic scales; use values on one scale.")
        if scales and scales != {'linear'}:
            metadata = dict(self.df.attrs.get(TARGET_SCALES_ATTR, {}))
            metadata[self.target_col] = next(iter(scales))
            self.df.attrs[TARGET_SCALES_ATTR] = metadata

    def _ensure_dtypes_after_load(self) -> None:
        """
        Ensure critical columns use stable pandas extension dtypes to avoid
        FutureWarning when assigning incompatible values.
        """
        # string columns
        for col in ["Canonical SMILES", "Processing Log", "Standardization Status"]:
            if col not in self.df.columns:
                self.df[col] = pd.Series(pd.NA, index=self.df.index, dtype="string")
            else:
                # Cast to pandas StringDtype regardless of current dtype
                self.df[col] = self.df[col].astype("string")

        # boolean column (nullable boolean dtype)
        if "Is Valid" not in self.df.columns:
            self.df["Is Valid"] = pd.Series(False, index=self.df.index, dtype="boolean")
        else:
            # Safe cast; non-boolean will be coerced to NA where needed
            try:
                self.df["Is Valid"] = self.df["Is Valid"].astype("boolean")
            except Exception:
                # Fallback: map common truthy/falsey to boolean, else NA
                self.df["Is Valid"] = (
                    self.df["Is Valid"]
                    .map(lambda v: True if v in [True, 1, "True", "true"] else (False if v in [False, 0, "False", "false"] else pd.NA))
                    .astype("boolean")
                )

    @check_data_loaded
    def preprocess(self, remove_salts: bool = True,
                   remove_solvents: bool = True,
                   remove_mixtures: bool = True,
                   remove_inorganic: bool = True,
                   reject_metal_species: bool = False,
                   element_policy: Optional[str] = None,
                   neutralize: bool = True,
                   reject_non_neutral: bool = False,
                   check_valid_atoms: bool = False,
                   remove_stereo: bool = False,
                   remove_isotopes: bool = True,
                   remove_hs: bool = True,
                   keep_largest_fragment: bool = False,
                   mixture_mode: Optional[str] = None,
                   hac_threshold: int = 3,
                   sanitize: bool = True,
                   reject_radical_species: bool = True,
                   progress_callback: Optional[Callable[[int, int], None]] = None,
                   n_jobs: int = 1,
                   chunksize: int = 100,
                   add_hs: bool = False) -> pd.DataFrame:
        """
        Execute the chemical processing pipeline.\
        :param remove_salts: Whether to remove salts.
        :param remove_solvents: Whether to remove solvent molecules.
        :param remove_mixtures: Legacy switch. False preserves mixtures; True uses reject/largest behavior.
        :param remove_inorganic: Whether to reject carbon-free and selected small inorganic structures.
        :param reject_metal_species: Whether to reject structures containing metal atoms.
        :param element_policy: Explicit mutually exclusive element policy: 'allow_all',
                               'reject_metals', or 'allowed_atoms'. When supplied, it takes
                               precedence over reject_metal_species and check_valid_atoms.
                               'allowed_atoms' rejects the whole molecule if any element is unlisted.
        :param neutralize: Whether to neutralize charges.
        :param reject_non_neutral: Only retain the molecules whose formal charge is zero.
        :param check_valid_atoms: Reject the whole molecule if any element is outside the allowed list.
        :param remove_stereo: Whether to remove stereochemistry.
        :param remove_isotopes: Whether to remove isotope information. Defaults to True.
        :param remove_hs: Whether to remove hydrogen atoms.
        :param add_hs: Add explicit hydrogen atoms after all chemical processing.
                       If remove_hs is also True, remove first and add last.
        :param keep_largest_fragment: Explicitly extract the largest substantial fragment instead of rejecting a mixture.
        :param mixture_mode: Explicit mixture policy: 'keep', 'reject', or 'largest'. When supplied,
                             it takes precedence over remove_mixtures and keep_largest_fragment.
                             All-identical repeated components collapse before any mixture policy.
        :param hac_threshold: Minimum heavy-atom threshold used only when extracting the largest fragment.
        :param sanitize: Whether to perform chemical sanitization.
        :param reject_radical_species: Molecules containing free radical atoms are directly rejected.
        :param progress_callback: Optional callback function for progressing.
        :param n_jobs: Number of parallel jobs to run.
                       1 means sequential execution (default).
                       -1 means use all available cores.
                       >1 means use specified number of cores.
        :param chunksize: Number of items to process per batch in multiprocessing.
        :return: Processed DataFrame with results.
        """
        if platform.system() == "Windows" and mp.current_process().name != 'MainProcess':
            msg = f"[WARNING] Process {os.getpid()} ignores 'preprocess' to prevent crash. Please use 'if __name__ == \"__main__\":' to fix memory issues."
            print(msg, file=sys.stderr, flush=True)
            return self.df

        mixture_aliases = {"preserve": "keep", "remove": "reject"}
        if mixture_mode is None:
            resolved_mixture_mode = (
                "keep" if not remove_mixtures
                else ("largest" if keep_largest_fragment else "reject")
            )
        else:
            resolved_mixture_mode = mixture_aliases.get(
                str(mixture_mode).strip().lower(),
                str(mixture_mode).strip().lower(),
            )
        if resolved_mixture_mode not in {"keep", "reject", "largest"}:
            raise ValueError("mixture_mode must be 'keep', 'reject', or 'largest'.")

        element_aliases = {
            "allow": "allow_all",
            "all": "allow_all",
            "metals": "reject_metals",
            "whitelist": "allowed_atoms",
        }
        if element_policy is None:
            resolved_element_policy = (
                "allowed_atoms" if check_valid_atoms
                else ("reject_metals" if reject_metal_species else "allow_all")
            )
        else:
            normalized_element_policy = str(element_policy).strip().lower()
            resolved_element_policy = element_aliases.get(
                normalized_element_policy,
                normalized_element_policy,
            )
        if resolved_element_policy not in {"allow_all", "reject_metals", "allowed_atoms"}:
            raise ValueError(
                "element_policy must be 'allow_all', 'reject_metals', or 'allowed_atoms'."
            )
        resolved_reject_metals = resolved_element_policy == "reject_metals"
        resolved_check_atoms = resolved_element_policy == "allowed_atoms"

        self._save_checkpoint()
        self._protect_source_columns(PREPROCESS_COLUMNS)
        self._ensure_dtypes_after_load()
        df_start = self.df.copy()

        # Capture configuration for the worker
        config = {
            'remove_salts': remove_salts,
            'remove_solvents': remove_solvents,
            'remove_mixtures': remove_mixtures,
            'mixture_mode': resolved_mixture_mode,
            'remove_inorganic': remove_inorganic,
            'reject_metal_species': resolved_reject_metals,
            'neutralize': neutralize,
            'reject_non_neutral': reject_non_neutral,
            'check_valid_atoms': resolved_check_atoms,
            'remove_stereo': remove_stereo,
            'remove_isotopes': remove_isotopes,
            'remove_hs': remove_hs,
            'add_hs': add_hs,
            'keep_largest_fragment': keep_largest_fragment,
            'hac_threshold': hac_threshold,
            'sanitize': sanitize,
            'reject_radical_species': reject_radical_species
        }

        smiles_list = self.df[self.smiles_col].tolist()
        total_rows = len(smiles_list)

        # Determine number of workers
        if n_jobs == -1:
            try:
                n_workers = mp.cpu_count()
            except:
                n_workers = 1
        else:
            n_workers = n_jobs

        results = []
        is_gui_mode = os.environ.get("DIPTOX_GUI_MODE") == "true"
        run_sequentially = False

        # Logic for Parallel vs Sequential
        if n_workers > 1:
            system_platform = platform.system()
            ctx = None

            try:
                # 1. Select optimal context
                if system_platform != "Windows":
                    # Linux/Mac: 'fork' is faster and doesn't require main guard
                    try:
                        ctx = mp.get_context('fork')
                    except ValueError:
                        ctx = mp.get_context('spawn')  # Fallback if fork unavailable
                else:
                    # Windows: Must use 'spawn'
                    ctx = mp.get_context('spawn')

                logger.info(f"Starting preprocessing with {n_workers} processes (OS: {system_platform})...")

                # Prepare arguments
                process_args = [(s, self.chem_processor, config) for s in smiles_list]

                # 2. Attempt execution
                with ctx.Pool(processes=n_workers) as pool:
                    iterator = pool.imap(
                        _worker_preprocess,
                        process_args,
                        chunksize=chunksize
                    )

                    for i, result in tqdm(enumerate(iterator), total=total_rows, desc=f"Processing (MP={n_workers})",
                                          disable=not self.interactive or bool(progress_callback)):
                        results.append(result)
                        if progress_callback and i % chunksize == 0:
                            progress_callback(i + 1, total_rows)

            except (RuntimeError, EOFError, OSError) as e:
                if mp.current_process().name == 'MainProcess':
                    err_msg = str(e).lower()
                    if "freeze_support" in err_msg or "pipe" in err_msg or "eof" in err_msg or system_platform == "Windows":
                        logger.error(
                            "Multiprocessing failed (likely due to missing main guard). Switching to sequential mode.")
                        if not is_gui_mode:
                            print(
                                "\n[System] Multiprocessing failed. Auto-fallback to sequential execution (n_jobs=1).\n",
                                file=sys.stderr)
                    else:
                        logger.error(f"Multiprocessing error: {e}")
                results = []
                run_sequentially = True
            except Exception as e:
                # Catch other unforeseen MP errors to prevent data loss
                if mp.current_process().name == 'MainProcess':
                    logger.error(f"Unexpected multiprocessing error: {e}. Fallback to sequential.")
                results = []
                run_sequentially = True
        else:
            run_sequentially = True

            # --- Sequential Fallback / Execution ---
        if run_sequentially:
            if n_workers > 1:
                logger.info("Running in sequential mode (Fallback)...")
            else:
                logger.info("Running in sequential mode...")

            for i, s in tqdm(enumerate(smiles_list), total=total_rows, desc="Processing",
                             disable=not self.interactive or bool(progress_callback)):
                result = _worker_preprocess((s, self.chem_processor, config))
                results.append(result)
                if progress_callback and i % 10 == 0:
                    progress_callback(i + 1, total_rows)

        # Keep malformed SDF records distinguishable from ordinary missing input.
        if 'Import Error' in self.df.columns:
            for index, import_error in enumerate(self.df['Import Error']):
                if pd.notna(import_error) and str(import_error).strip() and not results[index][0]:
                    results[index] = (False, f'Invalid SMILES: {import_error}', None)

        # Update DataFrame
        is_valid_list = [r[0] for r in results]
        logs_list = [r[1] for r in results]
        canon_smiles_list = [r[2] for r in results]
        status_list = [_standardization_status(r[0], r[1]) for r in results]
        original_metadata = [_molecule_metadata(smiles, sanitize=sanitize) for smiles in smiles_list]
        final_metadata = [_molecule_metadata(smiles) for smiles in canon_smiles_list]

        self.df['Is Valid'] = pd.Series(is_valid_list, index=self.df.index, dtype="boolean")
        self.df['Processing Log'] = pd.Series(logs_list, index=self.df.index, dtype="string")
        self.df['Canonical SMILES'] = pd.Series(canon_smiles_list, index=self.df.index, dtype="string")
        for column in self._structure_derived_columns:
            if column in self.df.columns:
                dtype = 'boolean' if column.startswith('Substructure_') else 'string'
                self.df[column] = pd.Series(pd.NA, index=self.df.index, dtype=dtype)
        self.df['Standardization Status'] = pd.Series(status_list, index=self.df.index, dtype="string")
        self.df['Original Canonical SMILES'] = pd.Series(
            [item.get('smiles') for item in original_metadata], index=self.df.index, dtype="string"
        )
        self.df['Original Fragment Count'] = pd.Series(
            [item.get('fragments') for item in original_metadata], index=self.df.index, dtype="Int64"
        )
        self.df['Final Fragment Count'] = pd.Series(
            [item.get('fragments') for item in final_metadata], index=self.df.index, dtype="Int64"
        )
        self.df['Original Structure Type'] = pd.Series(
            [
                "Multi-component" if item.get('fragments', 0) > 1 else "Single-component"
                if item else pd.NA
                for item in original_metadata
            ],
            index=self.df.index,
            dtype="string",
        )
        self.df['Final Structure Type'] = pd.Series(
            [
                "Multi-component" if item.get('fragments', 0) > 1 else "Single-component"
                if item else pd.NA
                for item in final_metadata
            ],
            index=self.df.index,
            dtype="string",
        )
        self.df['Is Multi-Component'] = pd.Series(
            [item.get('fragments', 0) > 1 if item else pd.NA for item in final_metadata],
            index=self.df.index,
            dtype="boolean",
        )
        self.df['Original Formal Charge'] = pd.Series(
            [item.get('charge') for item in original_metadata], index=self.df.index, dtype="Int64"
        )
        self.df['Final Formal Charge'] = pd.Series(
            [item.get('charge') for item in final_metadata], index=self.df.index, dtype="Int64"
        )
        self.df['Original Contains Metal'] = pd.Series(
            [bool(item.get('metals')) if item else pd.NA for item in original_metadata],
            index=self.df.index,
            dtype="boolean",
        )
        self.df['Final Contains Metal'] = pd.Series(
            [bool(item.get('metals')) if item else pd.NA for item in final_metadata],
            index=self.df.index,
            dtype="boolean",
        )
        self.df['Original Metal Elements'] = pd.Series(
            [",".join(item.get('metals', [])) if item else pd.NA for item in original_metadata],
            index=self.df.index,
            dtype="string",
        )
        self.df['Final Metal Elements'] = pd.Series(
            [",".join(item.get('metals', [])) if item else pd.NA for item in final_metadata],
            index=self.df.index,
            dtype="string",
        )
        self.df['Original Molecular Weight'] = pd.Series(
            [item.get('mw') for item in original_metadata], index=self.df.index, dtype="Float64"
        )
        self.df['Standardized Molecular Weight'] = pd.Series(
            [item.get('mw') for item in final_metadata], index=self.df.index, dtype="Float64"
        )
        self._invalidate_standardized_mw_targets()

        # Keep mapped source strings for provenance, but compare chemical
        # identities without atom-map annotations.
        original_identities = pd.Series(
            [item.get('identity') for item in original_metadata], index=self.df.index, dtype='string'
        )
        final_identities = pd.Series(
            [item.get('identity') for item in final_metadata], index=self.df.index, dtype='string'
        )
        self.df['Structure Changed'] = pd.Series(
            original_identities.ne(final_identities),
            index=self.df.index,
            dtype="boolean"
        )
        mw_delta = (self.df['Original Molecular Weight'] - self.df['Standardized Molecular Weight']).abs()
        self.df['Molecular Weight Changed'] = pd.Series(
            mw_delta.gt(1e-6), index=self.df.index, dtype="boolean"
        )

        self.df['Parent Mapping Count'] = pd.Series(pd.NA, index=self.df.index, dtype="Int64")
        self.df['Parent Mapping Conflict'] = pd.Series(pd.NA, index=self.df.index, dtype="boolean")
        self.df['Parent Source Structures'] = pd.Series(pd.NA, index=self.df.index, dtype="string")
        self.df['Parent Mapping Status'] = pd.Series("Excluded", index=self.df.index, dtype="string")
        valid_rows = self.df['Is Valid'].fillna(False) & self.df['Canonical SMILES'].notna()
        for _, group in self.df.loc[valid_rows].groupby('Canonical SMILES', sort=False):
            sources = list(dict.fromkeys(original_identities.loc[group.index].dropna().astype(str)))
            count = len(sources)
            indices = group.index
            self.df.loc[indices, 'Parent Mapping Count'] = count
            self.df.loc[indices, 'Parent Mapping Conflict'] = count > 1
            self.df.loc[indices, 'Parent Source Structures'] = " | ".join(sources)
            changed = bool(group['Structure Changed'].fillna(False).any())
            mapping_status = "Many-to-one" if count > 1 else ("Transformed" if changed else "Unchanged")
            self.df.loc[indices, 'Parent Mapping Status'] = mapping_status

        if progress_callback:
            progress_callback(total_rows, total_rows)
        self._preprocess_key = 1

        # Keep invalid rows in the working frame for inspection/reprocessing, and
        # expose an audit copy immediately instead of waiting for deduplication.
        excluded = self.df.loc[~self.df['Is Valid'].fillna(False)].copy()
        excluded['Source DataFrame Index'] = self._source_labels(excluded.index)
        excluded['Exclusion Reason'] = (
            excluded['Standardization Status'].fillna('Invalid').astype(str)
            + ': ' + excluded['Processing Log'].fillna('').astype(str)
        )
        excluded['Excluded At'] = 'Preprocessing'
        self._merge_excluded_records(excluded, refresh_frame=self.df, refresh_stage='Preprocessing')

        valid_count = self.df['Is Valid'].sum()
        invalid_count = total_rows - valid_count

        mode_str = "Multiprocessing" if (n_workers > 1 and not run_sequentially) else "Sequential"
        active_rules = []
        active_rules.append("Rej_Dummy_Atoms")
        if remove_salts: active_rules.append("Rm_Salts")
        if remove_solvents: active_rules.append("Rm_Solvents")
        mixture_rule = (
            f"Largest(HAC>{hac_threshold})"
            if resolved_mixture_mode == "largest"
            else resolved_mixture_mode.capitalize()
        )
        active_rules.append(f"Mixtures({mixture_rule})")
        if resolved_reject_metals: active_rules.append("Element_Policy(Reject_Metals)")
        if remove_inorganic: active_rules.append("Rm_Inorg")
        if neutralize: active_rules.append(f"Neutralize(Strict={reject_non_neutral})")
        if resolved_check_atoms: active_rules.append("Element_Policy(Allowed_Atoms)")
        if resolved_element_policy == "allow_all": active_rules.append("Element_Policy(Allow_All)")
        if remove_stereo: active_rules.append("Rm_Stereo")
        if remove_isotopes: active_rules.append("Rm_Iso")
        if remove_hs: active_rules.append("Rm_Hs")
        if add_hs: active_rules.append("Add_Hs")
        if sanitize: active_rules.append("Sanitize")
        if reject_radical_species: active_rules.append("Rej_Radicals")

        rules_str = ", ".join(active_rules) if active_rules else "None"

        details = f"Valid: {valid_count} | Invalid: {invalid_count} | Rules: [{rules_str}] | Mode: {mode_str}"
        self._record_step("Preprocessing", df_start, self.df, details)
        return self.df

    def _merge_excluded_records(self, excluded: pd.DataFrame,
                               refresh_frame: Optional[pd.DataFrame] = None,
                               refresh_stage: Optional[str] = None) -> None:
        """Refresh by private row identity, even for identical rows or duplicate labels."""
        previous = self.excluded_df
        if not previous.empty and refresh_frame is not None and not refresh_frame.empty:
            replaced = previous.index.isin(refresh_frame.index)
            if refresh_stage is not None and 'Excluded At' in previous:
                replaced &= previous['Excluded At'].eq(refresh_stage).fillna(False).to_numpy(dtype=bool)
            previous = previous.loc[~replaced]
        if previous.empty:
            self.excluded_df = excluded.copy()
        elif excluded.empty:
            self.excluded_df = previous.copy()
        else:
            self.excluded_df = pd.concat([previous, excluded])

    def _update_row(self, idx, is_valid: bool, comment: str, smiles: Optional[str]) -> None:
        """Update the result row for a given index."""
        self.df.at[idx, 'Is Valid'] = bool(is_valid)
        self.df.at[idx, 'Processing Log'] = "" if comment is None else str(comment)
        self.df.at[idx, 'Canonical SMILES'] = (pd.NA if smiles is None else str(smiles))

    @_run_on_main_process_only
    @check_data_loaded
    def standardize_units(self, standard_unit: Optional[str] = None,
                          conversion_rules: Optional[Dict[Tuple[str, str], str]] = None,
                          molecular_weight_source: Optional[str] = None,
                          molecular_weight_col: Optional[str] = None) -> None:
        """
        Orchestrates the standardization of units for the target column.
        :param standard_unit: The target unit to convert all values to.
        :param conversion_rules: A dictionary of conversion rules, e.g., {('mg/L', 'ug/L'): 'x * 1000'}.
        :param molecular_weight_source: Explicit molecular-weight basis for mass/molar conversion:
                                        'original' or 'standardized'. Required when a conversion
                                        uses MW, unless molecular_weight_col is provided.
        :param molecular_weight_col: Optional dataset column containing the reported molecular weight.
                                     When provided, it overrides molecular_weight_source.
        """
        df_start = self.df.copy()
        if molecular_weight_source not in {None, "original", "standardized"}:
            raise ValueError(
                "molecular_weight_source must be 'original' or 'standardized'; "
                "'auto' is no longer supported. Choose an explicit molecular-weight basis "
                "or provide molecular_weight_col."
            )
        if molecular_weight_source == "standardized" and not self._preprocess_key:
            raise ValueError(
                "molecular_weight_source='standardized' requires preprocessing first."
            )
        if not self.target_col or not self.unit_col:
            if not self.interactive:
                raise ValueError("Unit standardization requires both target_col and unit_col.")
            logger.info("Target column or unit column not specified, skipping unit standardization.")
            self._units_standardized = True
            return

        self._restore_target_scale_from_units()
        if get_target_scale(self.df, self.target_col) != 'linear':
            raise ValueError("Cannot convert a logarithmic target as a linear concentration. Undo the log transformation first.")
        molecular_weight_col = self._input_column_aliases.get(molecular_weight_col, molecular_weight_col)

        # Recompute from the last conversion's source columns on retry. DataFrame
        # attrs travel with checkpoints, so undo restores this mapping as well.
        source_target_col, source_unit_col = self.target_col, self.unit_col
        previous_conversion = self.df.attrs.get('_diptox_unit_conversion', {})
        conversions = self.df.attrs.get('_diptox_unit_conversions', {})
        if (
            previous_conversion.get('output_target') == self.target_col
            and previous_conversion.get('output_unit') == self.unit_col
            and previous_conversion.get('input_target') in self.df.columns
            and previous_conversion.get('input_unit') in self.df.columns
            and (standard_unit is None or standard_unit == previous_conversion.get('standard_unit'))
        ):
            source_target_col = previous_conversion['input_target']
            source_unit_col = previous_conversion['input_unit']
            if standard_unit is None:
                standard_unit = previous_conversion.get('standard_unit')

        # A repeated request for the current unit can change the MW basis even
        # after intermediate unit scaling. Find the actual MW conversion's input.
        requested_basis = f'column:{molecular_weight_col}' if molecular_weight_col else molecular_weight_source
        current_units = self.df[self.unit_col].dropna().unique()
        same_unit_request = len(current_units) == 1 and (standard_unit is None or standard_unit == current_units[0])
        if requested_basis and same_unit_request:
            cursor, visited = self.target_col, set()
            while cursor in conversions and cursor not in visited:
                visited.add(cursor)
                conversion = conversions[cursor]
                if conversion.get('requires_mw'):
                    if (requested_basis != conversion.get('mw_basis')
                            or self.target_col in self.df.attrs.get('_diptox_stale_targets', [])):
                        source_target_col = conversion['input_target']
                        source_unit_col = conversion['input_unit']
                    break
                cursor = conversion['input_target']
            else:
                prior_basis = self.df.attrs.get('_diptox_target_mw_bases', {}).get(self.target_col)
                if prior_basis and prior_basis != requested_basis:
                    raise ValueError("Changing the molecular-weight basis after aggregation requires restoring the original concentrations first.")

        if source_target_col in self.df.attrs.get('_diptox_stale_targets', []):
            raise ValueError("The target depends on an outdated standardized molecular weight. Restore the original concentrations before converting again.")

        unique_units = [u for u in self.df[source_unit_col].dropna().unique() if u]
        current_unit = unique_units[0] if unique_units else None
        final_standard_unit = standard_unit or (current_unit if len(unique_units) <= 1 else None)
        is_gui_mode = os.environ.get("DIPTOX_GUI_MODE") == "true"

        if not final_standard_unit:
            if not self.interactive or is_gui_mode:
                raise ValueError(
                    "A standard_unit must be provided when source units are missing or multiple units exist. "
                    f"Available source units: {unique_units}"
                )
            final_standard_unit = self._select_standard_unit_interactively(unique_units)
            if not final_standard_unit:
                return

        unit_processor = UnitProcessor(rules=conversion_rules)

        if not self.interactive:
            missing_rules = [
                (unit, final_standard_unit) for unit in unique_units
                if unit != final_standard_unit and not unit_processor.get_rule(unit, final_standard_unit)
            ]
            if missing_rules:
                pairs = ", ".join(f"{source!r} -> {target!r}" for source, target in missing_rules)
                raise ValueError(
                    f"Missing conversion rules: {pairs}. Provide conversion_rules with a formula for each pair."
                )

        rule_provider = None
        if self.interactive and not is_gui_mode:
            prompt_tracker = {'first_time': True}
            rule_provider = lambda from_unit, to_unit: self._get_rule_from_user(from_unit, to_unit, prompt_tracker)

        smiles_col = 'Canonical SMILES' if molecular_weight_source == "standardized" else None
        mw_col = molecular_weight_col
        if not mw_col and molecular_weight_source == "original":
            if self._preprocess_key and 'Original Molecular Weight' in self.df.columns:
                mw_col = 'Original Molecular Weight'
            else:
                smiles_col = self.smiles_col

        requires_mw = any(
            unit != final_standard_unit
            and 'mw' in (unit_processor.get_rule(unit, final_standard_unit) or '')
            for unit in unique_units
        )
        if requires_mw and not mw_col and molecular_weight_source is None:
            raise ValueError(
                "This conversion requires an explicit molecular-weight basis. Set "
                "molecular_weight_source='original' or 'standardized', "
                "or provide molecular_weight_col."
            )

        self._save_checkpoint()
        while True:
            renamed = self._protect_source_columns({
                source_target_col + ' (Standardized)', source_unit_col + ' (Standardized)',
                'Unit Conversion Status', 'Molecular Weight Used', 'Molecular Weight Source',
            })
            if not renamed:
                break
            source_target_col = renamed.get(source_target_col, source_target_col)
            source_unit_col = renamed.get(source_unit_col, source_unit_col)
            if mw_col:
                mw_col = renamed.get(mw_col, mw_col)
            if smiles_col:
                smiles_col = renamed.get(smiles_col, smiles_col)
        if molecular_weight_col:
            requested_basis = f'column:{mw_col}'
        try:
            converted_df, new_target_col, new_unit_col = unit_processor.standardize(
                df=self.df,
                target_col=source_target_col,
                unit_col=source_unit_col,
                standard_unit=final_standard_unit,
                smiles_col=smiles_col,
                molecular_weight_col=mw_col,
                rule_provider_callback=rule_provider
            )
            converted_df.attrs['_diptox_unit_conversion'] = {
                'input_target': source_target_col,
                'input_unit': source_unit_col,
                'output_target': new_target_col,
                'output_unit': new_unit_col,
                'standard_unit': final_standard_unit,
                'requires_mw': requires_mw,
                'mw_basis': requested_basis,
            }
            conversion_history = deepcopy(self.df.attrs.get('_diptox_unit_conversions', {}))
            conversion_history[new_target_col] = dict(converted_df.attrs['_diptox_unit_conversion'])
            converted_df.attrs['_diptox_unit_conversions'] = conversion_history
            bases = dict(self.df.attrs.get('_diptox_target_mw_bases', {}))
            basis = requested_basis if requires_mw else bases.get(source_target_col)
            if basis:
                bases[new_target_col] = basis
            converted_df.attrs['_diptox_target_mw_bases'] = bases
            dependencies = deepcopy(self.df.attrs.get('_diptox_standardized_mw_dependencies', {}))
            inherited = dict(dependencies.get(source_target_col, {}))
            if requires_mw and ((molecular_weight_source == 'standardized' and not mw_col)
                                or (mw_col == 'Standardized Molecular Weight' and self._preprocess_key)):
                for index in converted_df.index:
                    used = converted_df.at[index, 'Molecular Weight Used']
                    if pd.notna(used):
                        inherited[index] = float(used)
            if inherited:
                dependencies[new_target_col] = inherited
            else:
                dependencies.pop(new_target_col, None)
            converted_df.attrs['_diptox_standardized_mw_dependencies'] = dependencies
            converted_df.attrs['_diptox_stale_targets'] = [
                column for column in self.df.attrs.get('_diptox_stale_targets', [])
                if column != new_target_col
            ]
            self.df = converted_df
            self.target_col = new_target_col
            self.unit_col = new_unit_col
            self._units_standardized = True
            mw_basis = f"column:{mw_col}" if mw_col else molecular_weight_source or "not required"
            self._record_step(
                "Unit Standardization", df_start, self.df,
                f"Target: {final_standard_unit} | MW basis: {mw_basis}"
            )
        except ValueError as e:
            logger.error(f"Unit standardization failed: {e}")
            raise

    @staticmethod
    def _select_standard_unit_interactively(unique_units: List[str]) -> Optional[str]:
        """Handles the interactive command-line prompt for selecting a standard unit."""
        print("Multiple units found. Please select a standard unit to convert to:")
        for i, unit in enumerate(unique_units):
            print(f"{i + 1}: {unit}")

        while True:
            try:
                choice_input = input(f"Enter the number of the standard unit (1-{len(unique_units)}): ").strip()
                choice = int(choice_input) - 1
                if 0 <= choice < len(unique_units):
                    selected_unit = unique_units[choice]
                    logger.info(f"Standard unit set to '{selected_unit}'.")
                    return selected_unit
                else:
                    print("Invalid number. Please try again.")
            except (ValueError, IndexError):
                print("Invalid input. Please enter a number.")
            except (EOFError, KeyboardInterrupt):
                logger.warning("Unit selection cancelled by user.")
                return None

    @staticmethod
    def _get_rule_from_user(from_unit: str, to_unit: str, tracker: Dict) -> Optional[str]:
        """Callback function to get conversion rules from the command-line user."""
        print("-" * 30)
        print(f"No conversion rule found for '{from_unit}' -> '{to_unit}'.")

        if tracker['first_time']:
            print("Please provide a conversion formula (use 'x' for the value).")
            print("Example (mg/L to ug/L): x * 1000")
            print("Example (log(mol/L) to mol/L): 10**(-x)")
            tracker['first_time'] = False

        try:
            prompt = "Enter formula (press Enter to skip; data with this unit will be removed): "
            formula_input = input(prompt).strip()
            return formula_input if formula_input else None
        except (EOFError, KeyboardInterrupt):
            logger.warning(f"Rule input for '{from_unit}' cancelled.")
            return None

    @_run_on_main_process_only
    def config_deduplicator(self, condition_cols: Optional[List[str]] = None,
                            data_type: str = "continuous",
                            method: str = "auto",
                            priority: Optional[List[str]] = None,
                            p_threshold: float = 0.05,
                            custom_method: Optional[Callable] = None,
                            standard_unit: Optional[str] = None,
                            conversion_rules: Optional[Dict[Tuple[str, str], str]] = None,
                            molecular_weight_source: Optional[str] = None,
                            molecular_weight_col: Optional[str] = None,
                            log_transform: Union[bool, str] = "None",
                            dropna_conditions: bool = False,
                            aggregation: str = "mean") -> None:
        """
        Configure the deduplicator device
        :param condition_cols: Experimental context columns that define distinct modeling records
                               (e.g. temperature or species). Do not include removed salt/metal
                               provenance unless it is also an input feature of the final model.
        :param data_type: data type - discrete/continuous
        :param method: Continuous outlier filtering (auto, 3sigma, IQR), or discrete
                       selection (vote, priority). Groups of <=3 skip built-in filtering.
        :param aggregation: Continuous aggregation (mean, max, min) after transformation
                            and filtering. Extrema ties retain the first matching record.
        :param priority: Ordered preferred values, used only by discrete method='priority'.
                         Other methods (including vote) ignore it. No match falls back to voting.
        :param p_threshold: Threshold of normal distribution
        :param custom_method: Custom method of data deduplication
        :param standard_unit: The target unit to standardize to before deduplication.
        :param conversion_rules: A dictionary of rules for unit conversion.
        :param molecular_weight_source: Molecular-weight basis used by an implicit unit conversion.
        :param molecular_weight_col: Optional explicit molecular-weight column.
        :param log_transform: If True, applies a -log10 transformation to the target column.
        :param dropna_conditions: If True, drops rows with missing condition values. If False, groups them.
        """
        if molecular_weight_source not in {None, "original", "standardized"}:
            raise ValueError(
                "molecular_weight_source must be 'original' or 'standardized'; "
                "'auto' is no longer supported."
            )
        self._dedup_unit_settings = (
            {
                'standard_unit': standard_unit,
                'conversion_rules': conversion_rules,
                'molecular_weight_source': molecular_weight_source,
                'molecular_weight_col': molecular_weight_col,
            }
            if standard_unit or conversion_rules else None
        )
        if condition_cols:
            condition_cols = [self._input_column_aliases.get(column, column) for column in condition_cols]

        smiles_col = 'Canonical SMILES' if self._preprocess_key else self.smiles_col

        self.deduplicator = DataDeduplicator(
            smiles_col=smiles_col, target_col=self.target_col, condition_cols=condition_cols,
            data_type=data_type, method=method, p_threshold=p_threshold, priority=priority,
            custom_method=custom_method, log_transform=log_transform, dropna_conditions=dropna_conditions,
            aggregation=aggregation
        )
        self._current_dedup_config = {
            'method': method,
            'aggregation': aggregation,
            'data_type': data_type,
            'condition_cols': condition_cols,
            'priority': self.deduplicator.priority_list,
            'log_transform': log_transform,
            'dropna_conditions': dropna_conditions
        }

    @_run_on_main_process_only
    @check_data_loaded
    def dataset_deduplicate(self, progress_callback: Optional[Callable] = None) -> None:
        """Execution deduplicator removal"""
        if not self.deduplicator:
            raise ValueError("Deduplicator not configured. Call config_deduplicator first.")
        self._restore_target_scale_from_units()
        self._save_checkpoint()
        df_start = self.df.copy()

        # Resolve the structure key at execution time. This guarantees that a
        # deduplicator configured before preprocessing still uses the final
        # representation produced by the user's selected preprocessing policy.
        dedup_structure_col = 'Canonical SMILES' if self._preprocess_key else self.smiles_col
        if not dedup_structure_col or dedup_structure_col not in self.df.columns:
            raise ValueError("No structure column is available for deduplication.")
        self.deduplicator.smiles_col = dedup_structure_col

        if self._dedup_unit_settings:
            logger.info("Implicitly running unit standardization as part of deduplication.")
            self.standardize_units(
                standard_unit=self._dedup_unit_settings.get('standard_unit'),
                conversion_rules=self._dedup_unit_settings.get('conversion_rules'),
                molecular_weight_source=self._dedup_unit_settings.get('molecular_weight_source'),
                molecular_weight_col=self._dedup_unit_settings.get('molecular_weight_col')
            )
            self._dedup_unit_settings = None
            df_start = self.df.copy()

        if self.deduplicator.data_type != 'smiles' and self.target_col in self.df.attrs.get('_diptox_stale_targets', []):
            raise ValueError("The target uses an outdated standardized molecular weight. Reconvert from the original concentrations before deduplication.")

        needs_standardization = False
        if self.target_col and self.unit_col and self.unit_col in self.df.columns:
            unique_count = self.df[self.unit_col].dropna().nunique()
            needs_standardization = unique_count > 1

        if needs_standardization and not self._units_standardized:
            raise ValueError(
                "Unit standardization is required but has not been performed. Please go to the 'Unit Standardization' step first.")

        input_target = self.target_col
        input_scale = get_target_scale(self.df, input_target)
        while True:
            outputs = {'Deduplication Strategy', 'Deduplication Record Count', 'Deduplication Source Rows',
                       'Deduplication Input Values', 'Deduplication Distinct Value Count', 'Deduplication Value Range'}
            if self.target_col:
                outputs.add(self.target_col + '_new')
            if not self._protect_source_columns(outputs):
                break
        input_target = self.target_col
        self.deduplicator.target_col = input_target
        dependencies = deepcopy(self.df.attrs.get('_diptox_standardized_mw_dependencies', {}))

        self.df = self.deduplicator.deduplicate(self.df, progress_callback=progress_callback)
        for index, members in self.deduplicator.source_indices.items():
            original_members = tuple(dict.fromkeys(
                source for member in members for source in self._row_lineage.get(member, (member,))
            ))
            self._row_lineage[index] = original_members
            self.df.at[index, 'Deduplication Source Rows'] = ';'.join(
                str(label) for label in self._source_labels(original_members)
            )
            self.df.at[index, 'Deduplication Record Count'] = len(original_members)
        exclusion_reasons = self.deduplicator.exclusion_reasons
        excluded_indices = [index for index in df_start.index if index in exclusion_reasons]
        if excluded_indices:
            excluded = df_start.loc[excluded_indices].copy()
            reasons = []
            for index, row in excluded.iterrows():
                reason = exclusion_reasons[index]
                std_status = row.get('Standardization Status')
                if pd.notna(std_status) and std_status != "Retained":
                    processing_log = row.get('Processing Log')
                    reason = str(std_status)
                    if pd.notna(processing_log) and str(processing_log).strip():
                        reason += f": {processing_log}"
                else:
                    unit_status = row.get('Unit Conversion Status')
                    successful_unit_statuses = {"Converted", "No conversion needed"}
                    if pd.notna(unit_status) and unit_status not in successful_unit_statuses:
                        reason = f"Unit conversion: {unit_status}"
                reasons.append(reason)
            excluded['Deduplication Exclusion Reason'] = reasons
            excluded['Exclusion Reason'] = reasons
            excluded['Excluded At'] = "Deduplication"
            excluded['Source DataFrame Index'] = self._source_labels(excluded.index)
            self._merge_excluded_records(excluded, refresh_frame=excluded)
        if self.target_col and self.deduplicator.data_type != 'smiles' and self.target_col + '_new' in self.df:
            self.target_col = self.target_col + "_new"
            input_weights = dependencies.get(input_target, {})
            propagated = {}
            for index, members in self.deduplicator.source_indices.items():
                values = [input_weights[member] for member in members if member in input_weights]
                if values:
                    propagated[index] = values[0]
            if propagated:
                dependencies[self.target_col] = propagated
            self.df.attrs['_diptox_standardized_mw_dependencies'] = dependencies
            bases = dict(self.df.attrs.get('_diptox_target_mw_bases', {}))
            if input_target in bases:
                bases[self.target_col] = bases[input_target]
            self.df.attrs['_diptox_target_mw_bases'] = bases
            output_scale = get_target_scale(self.df, self.target_col)
            if output_scale != input_scale and self.unit_col and self.unit_col in self.df:
                output_unit = self.unit_col + '_new'
                self._protect_source_columns({output_unit})
                source_units = self.df[self.unit_col].astype('string')
                self.df[output_unit] = output_scale + '(' + source_units + ')'
                self.unit_col = output_unit

        method_name = self.deduplicator.method
        if self.deduplicator.log_transform:
            method_name += " (Log10 Transformed)"

        cfg = getattr(self, '_current_dedup_config', {})
        method_name = cfg.get('method', 'unknown')
        if cfg.get('data_type') == 'continuous':
            method_name += f" -> {cfg.get('aggregation', 'mean')}"
        trans = cfg.get('log_transform', "None")
        if trans != "None":
            method_name += f" ({trans} Transformed)"
        conds_list = cfg.get('condition_cols')
        conds = f"Conds: {conds_list}" if conds_list else "No Conds"
        dropna_str = "DropNA" if cfg.get('dropna_conditions', False) else "KeepNA"
        priority_list = cfg.get('priority')
        priority_str = f" | Priority: {priority_list}" if priority_list else ""
        excluded_count = len(self.excluded_df)
        details = (
            f"Method: {method_name} ({cfg.get('data_type', 'unknown')}) | "
            f"Structure key: {dedup_structure_col} | {conds} | {dropna_str}{priority_str} | "
            f"Excluded: {excluded_count}"
        )
        self._record_step("Deduplication", df_start, self.df, details)

    @_run_on_main_process_only
    @check_data_loaded
    def substructure_search(self, query_pattern: Union[str, List[str]],
                            is_smarts: bool = False) -> None:
        """
        Integrated search interface
        :param query_pattern: Molecular substructure
        :param is_smarts: Search mode (SMILES/SMARTS)
        """
        self._save_checkpoint()
        df_start = self.df.copy()
        query_pattern_list = [query_pattern] if isinstance(query_pattern, str) else query_pattern
        self._protect_source_columns({f'Substructure_{pattern}' for pattern in query_pattern_list})
        searcher = SubstructureSearcher(
            df=self.df,
            smiles_col='Canonical SMILES' if self._preprocess_key else self.smiles_col,
        )
        for query_pattern in query_pattern_list:
            results = searcher.search(query_pattern, is_smarts)
            col_name = f'Substructure_{query_pattern}'
            self.df[col_name] = pd.Series(False, index=self.df.index, dtype="boolean")
            self._structure_derived_columns.add(col_name)

            for idx, _ in results['matches']:
                self.df.at[idx, col_name] = True
            matches = self.df[col_name].sum()
            self._record_step("Substructure Search", df_start, self.df,
                              f"Pattern: {query_pattern} (SMARTS: {is_smarts}) | Matches Found: {matches}")

    @_run_on_main_process_only
    def config_web_request(self, sources: Union[str, List[str]] = 'pubchem',
                           interval: int = 0.3,
                           retries: int = 3,
                           delay: int = 30,
                           max_workers: int = 1,
                           batch_limit: int = 1500,
                           rest_duration: int = 300,
                           chemspider_api_key: Optional[str] = None,
                           comptox_api_key: Optional[str] = None,
                           cas_api_key: Optional[str] = None,
                           force_api_mode: bool = False,
                           status_callback: Optional[Callable[[str], None]] = None) -> None:
        """
        Initializes the WebService class
        :param sources: Data source interface
        :param interval: Time interval in seconds
        :param retries: Number of retry attempts on failure.
        :param delay: Delay between retries (in seconds).
        :param max_workers: Maximum number of concurrent requests.
        :param batch_limit: Number of requests before taking a break.
        :param rest_duration: Duration of the break in seconds.
        :param chemspider_api_key: Chemspider API key.
        :param comptox_api_key: Comptox API key.
        :param cas_api_key: CAS API key.
        :param force_api_mode: Force API mode.
        :param status_callback: Callback function for status checking.
        """
        self.web_source = sources
        self.web_service = WebService(sources=sources, interval=interval, retries=retries, delay=delay,
                                      max_workers=max_workers, batch_limit=batch_limit, rest_duration=rest_duration,
                                      chemspider_api_key=chemspider_api_key, comptox_api_key=comptox_api_key,
                                      cas_api_key=cas_api_key, force_api_mode=force_api_mode,
                                      status_callback=status_callback)

    @_run_on_main_process_only
    @check_data_loaded
    def web_request(self, send: Union[str, List[str]], request: Union[str, List[str]],
                    progress_callback: Optional[Callable] = None) -> None:
        """
        Add CAS numbers for valid molecules.
        :param send: What is used to request additional data? (smiles/cas)
        :param request: What identifier is requested? (smiles/cas/iupac)
        """
        if self.web_service is None:
            raise ValueError("The WebService has not been configured. Please call the config_web_request first.")
        self._save_checkpoint()
        df_start = self.df.copy()
        send_by_list = [send] if isinstance(send, str) else send
        send_ordered = list(dict.fromkeys([prop.strip().lower() for prop in send_by_list]))
        request_list = [request] if isinstance(request, str) else request
        request_set = {prop.strip().lower() for prop in request_list}

        VALID_PROPERTIES = {'smiles', 'cas', 'iupac', 'mw', 'name'}
        invalid_props = request_set - VALID_PROPERTIES
        if invalid_props:
            raise ValueError(f"Invalid request properties: {list(invalid_props)}")

        smiles_col = 'Canonical SMILES' if self._preprocess_key else self.smiles_col
        col_map = {'cas': self.cas_col, 'name': self.name_col, 'smiles': smiles_col}

        for prop in request_set:
            col = f'{prop}_from_web'
            if col not in self.df.columns:
                self.df[col] = pd.Series(pd.NA, index=self.df.index, dtype="string")
            else:
                try:
                    self.df[col] = self.df[col].astype("string")
                except Exception:
                    self.df[col] = pd.Series(pd.NA, index=self.df.index, dtype="string")

        for meta_col in ['Query_Status', 'Data_Source', 'Query_Method']:
            if meta_col not in self.df.columns:
                self.df[meta_col] = pd.Series(pd.NA, index=self.df.index, dtype="string")
            else:
                try:
                    self.df[meta_col] = self.df[meta_col].astype("string")
                except Exception:
                    self.df[meta_col] = pd.Series(pd.NA, index=self.df.index, dtype="string")

        self.df['Query_Status'] = "Pending"

        pending_indices = list(self.df.index)
        total_queries = len(pending_indices)

        try:
            total_steps = len(send_ordered)
            for step_idx, id_type in enumerate(send_ordered):
                if not pending_indices:
                    break

                input_col_name = col_map.get(id_type)
                if not input_col_name or input_col_name not in self.df.columns:
                    logger.warning(f"The specified column '{input_col_name}' (for '{id_type}') does not exist, so it will be skipped.")
                    continue

                identifiers_to_query = self.df.loc[pending_indices, input_col_name].tolist()

                if id_type == 'cas':
                    identifiers_to_query = [self.web_service._validate_and_clean_cas(cas) for cas in
                                            identifiers_to_query]

                results = self.web_service.get_properties_batch(
                    identifiers_to_query,
                    request_set,
                    id_type,
                    progress_callback=progress_callback
                )

                processed_indices = []
                for i, original_index in enumerate(pending_indices):
                    res = results[i]
                    if any(res.get(prop) for prop in request_set):
                        for prop in request_set:
                            col = f'{prop}_from_web'
                            val = res.get(prop)

                            if val is not None and pd.isna(self.df.at[original_index, col]):
                                self.df.at[original_index, col] = str(val)

                        self.df.at[original_index, 'Data_Source'] = res.get('Data_Source')
                        self.df.at[original_index, 'Query_Method'] = id_type
                        self.df.at[original_index, 'Query_Status'] = 'Success'
                        processed_indices.append(original_index)
                    else:
                        self.df.at[original_index, 'Data_Source'] = res.get('Data_Source')

                pending_indices = [idx for idx in pending_indices if idx not in processed_indices]

            if pending_indices:
                self.df.loc[pending_indices, 'Query_Status'] = 'Failed'

            logger.info("Web request processing is complete.")

            if 'smiles' in request_set:
                if not self._preprocess_key:
                    self.smiles_col = 'smiles_from_web'
            if 'cas' in request_set:
                self.cas_col = 'cas_from_web'
            if 'name' in request_set:
                self.name_col = 'name_from_web'

            success_count = len(self.df[self.df['Query_Status'] == 'Success'])
            props_str = ", ".join(request_set)
            send_str = " -> ".join(send_ordered)
            sources_str = ", ".join(self.web_source) if isinstance(self.web_source, list) else self.web_source

            details = f"Queried: {props_str} (via {send_str}) | Sources: [{sources_str}] | Success: {success_count}/{total_queries}"
            self._record_step("Web Request", df_start, self.df, details)

        except requests.exceptions.RequestException as e:
            logger.error(f"An error occurred on the network, and the web request was interrupted: {e}")
            self.df.loc[pending_indices, 'Query_Status'] = 'Network Error'
            return
        except Exception as e:
            logger.error(f"An unknown error occurred during the web request: {e}")
            self.df.loc[pending_indices, 'Query_Status'] = 'Unknown Error'
            raise

    @_run_on_main_process_only
    @check_data_loaded
    def calculate_inchi(self) -> None:
        """
        Calculate InChI strings locally using RDKit based on the current SMILES column.
        No web request required.
        """
        inchi_col = 'InChI'

        self._save_checkpoint()
        self._protect_source_columns({inchi_col})
        smiles_col = 'Canonical SMILES' if self._preprocess_key else self.smiles_col
        self.df[inchi_col] = pd.Series(pd.NA, index=self.df.index, dtype="string")
        self._structure_derived_columns.add(inchi_col)

        logger.info("Calculating InChI from SMILES using RDKit...")

        count = 0
        for idx, row in tqdm(self.df.iterrows(), total=len(self.df), desc="InChI Calc",
                             disable=not self.interactive):
            smiles = row[smiles_col]
            if not isinstance(smiles, str) or not smiles.strip():
                continue

            mol = self.chem_processor.smiles_to_mol(smiles, sanitize=True)
            inchi = self.chem_processor.mol_to_inchi(mol)

            if inchi:
                self.df.at[idx, inchi_col] = inchi
                count += 1

        logger.info(f"InChI calculation complete. Generated {count} InChI strings.")
        self._record_step("InChI Calculation", None, self.df, f"Calculated {count} InChIs")

    @_run_on_main_process_only
    @check_data_loaded
    def filter_by_atom_count(self,
                             min_heavy_atoms: Optional[int] = None,
                             max_heavy_atoms: Optional[int] = None,
                             min_total_atoms: Optional[int] = None,
                             max_total_atoms: Optional[int] = None) -> None:
        """
        Filter molecules based on heavy or total atom counts.
        :param min_heavy_atoms: Minimum number of heavy atoms (inclusive).
        :param max_heavy_atoms: Maximum number of heavy atoms (inclusive).
        :param min_total_atoms: Minimum number of total atoms (inclusive).
        :param max_total_atoms: Maximum number of total atoms (inclusive).
        """
        self._save_checkpoint()
        df_start = self.df.copy()
        if all(arg is None for arg in [min_heavy_atoms, max_heavy_atoms, min_total_atoms, max_total_atoms]):
            logger.warning("No filter criteria provided for filter_by_atom_count. No action taken.")
            return

        initial_count = len(self.df)
        smiles_col = 'Canonical SMILES' if self._preprocess_key else self.smiles_col

        def is_valid_by_atom_count(row):
            s = row[smiles_col]
            if not isinstance(s, str) or not s.strip():
                return False

            if self._preprocess_key and 'Is Valid' in row.index and not pd.isna(row['Is Valid']):
                if not row['Is Valid']:
                    return False

            mol = self.chem_processor.smiles_to_mol(s, sanitize=True)

            return self.chem_processor.validate_atom_count(
                mol,
                min_heavy_atoms,
                max_heavy_atoms,
                min_total_atoms,
                max_total_atoms
            )

        mask = self.df.apply(is_valid_by_atom_count, axis=1)
        newly_excluded = ~mask
        if self._preprocess_key:
            newly_excluded &= self.df['Is Valid'].fillna(False)
        excluded = self.df.loc[newly_excluded].copy()
        if not excluded.empty:
            excluded['Source DataFrame Index'] = self._source_labels(excluded.index)
            excluded['Excluded At'] = 'Atom count filter'
            excluded['Exclusion Reason'] = (
                f'Atom count outside requested range: heavy={min_heavy_atoms}-{max_heavy_atoms}, '
                f'total={min_total_atoms}-{max_total_atoms}'
            )
            self._merge_excluded_records(excluded, refresh_frame=excluded)
        self.df = self.df.loc[mask].copy()
        final_count = len(self.df)
        logger.info(
            f"Filtered by atom count. Initial: {initial_count}, Final: {final_count}, Removed: {initial_count - final_count}")
        details = f"Heavy: {min_heavy_atoms}-{max_heavy_atoms} | Total: {min_total_atoms}-{max_total_atoms} | Removed: {initial_count - final_count}"
        self._record_step("Filter Atom Count", df_start, self.df, details)

    @_run_on_main_process_only
    @check_data_loaded
    def save_results(self, output_path: str, columns: Optional[List[str]] = None) -> Optional[str]:
        """
        Save the processed results to a file.
        :param output_path: The output path where the results will be saved.
        :param columns: The columns to save (default saves all columns).
        """
        save_cols = columns if columns else self.df.columns.tolist()
        return self.data_handler.save_data(
            self.df, output_path, save_cols,
            'Canonical SMILES' if self._preprocess_key else self.smiles_col,
            self.id_col, interactive=self.interactive,
            source_smiles_col=self.smiles_col,
        )

    @_run_on_main_process_only
    def get_excluded_records(self) -> pd.DataFrame:
        """Return source rows that could not produce a modeling record during deduplication."""
        return self.excluded_df.copy()

    @_run_on_main_process_only
    def save_excluded_results(self, output_path: str, columns: Optional[List[str]] = None) -> Optional[str]:
        """Save deduplication exclusions with record-level reasons for audit."""
        if self.excluded_df.empty:
            raise ValueError("No excluded deduplication records are available.")
        if output_path.lower().endswith(('.sdf', '.smi')):
            raise ValueError("Excluded records must be exported as CSV, TXT, XLS, or XLSX.")
        save_cols = columns if columns else self.excluded_df.columns.tolist()
        missing_cols = [column for column in save_cols if column not in self.excluded_df.columns]
        if missing_cols:
            raise KeyError(f"Columns not found in excluded records: {missing_cols}")
        smiles_col = 'Canonical SMILES' if self._preprocess_key else self.smiles_col
        return self.data_handler.save_data(
            self.excluded_df,
            output_path,
            save_cols,
            smiles_col,
            self.id_col,
            interactive=self.interactive,
            source_smiles_col=self.smiles_col,
        )

    # Chemical rule management interface
    @_run_on_main_process_only
    def add_neutralization_rule(self, reactant: str, product: str) -> None:
        """
        Add a new neutralization rule to the list, ensuring the rule is valid and there are no conflicts.
        :param reactant: SMARTS string for the reactant.
        :param product: SMILES string for the product.
        """
        return self.chem_processor.add_neutralization_rule(reactant, product)

    @_run_on_main_process_only
    def remove_neutralization_rule(self, reactant: str) -> None:
        """
        Remove a charge neutralization rule.
        :param reactant: SMARTS string for the reactant.
        """
        return self.chem_processor.remove_neutralization_rule(reactant)

    @_run_on_main_process_only
    def manage_atom_rules(self, atoms: Union[str, List[str]], add: bool = True) -> List[str]:
        """Manage atom validation rules."""
        atom_list = [atoms] if isinstance(atoms, str) else atoms
        failed = []
        for atom in atom_list:
            if add:
                success = self.chem_processor.add_effective_atom(atom)
            else:
                success = self.chem_processor.delete_effective_atom(atom)
            if not success:
                failed.append(atom)
        return failed

    @_run_on_main_process_only
    def manage_default_salt(self, salts: Union[str, List[str]], add: bool = True) -> List[str]:
        """Manage salt validation rules."""
        salt_list = [salts] if isinstance(salts, str) else salts
        failed = []
        for salt in salt_list:
            if add:
                success = self.chem_processor.add_default_salt(salt)
            else:
                success = self.chem_processor.remove_default_salt(salt)
            if not success:
                failed.append(salt)
        return failed

    @_run_on_main_process_only
    def manage_default_solvent(self, solvents: Union[str, List[str]], add: bool = True) -> List[str]:
        """Manage solvent validation rules."""
        solvent_list = [solvents] if isinstance(solvents, str) else solvents
        failed = []
        for solvent in solvent_list:
            if add:
                success = self.chem_processor.add_default_solvents(solvent)
            else:
                success = self.chem_processor.remove_default_solvents(solvent)
            if not success:
                failed.append(solvent)
        return failed

    @_run_on_main_process_only
    def display_processing_rules(self) -> None:
        """
        Displays the current chemical processing rules being used,
        including valid atoms, neutralization rules, salts, and solvents.
        """
        self.chem_processor.display_current_rules()
