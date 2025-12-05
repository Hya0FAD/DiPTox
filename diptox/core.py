# diptox/core.py
import os
import pandas as pd
from typing import Optional, List, Union, Tuple, Callable, Dict
from functools import partial, wraps
from tqdm import tqdm
from rdkit import Chem
from datetime import datetime
import requests
from .chem_processor import ChemistryProcessor
from .web_request import WebService
from .data_io import DataHandler
from .data_deduplicator import DataDeduplicator
from .substructure_search import SubstructureSearcher
from .unit_processor import UnitProcessor
from .logger import log_manager
logger = log_manager.get_logger(__name__)
from diptox import user_reg


def check_data_loaded(func):
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        if self.df is None:
            raise ValueError("No data loaded. Please load data first.")
        return func(self, *args, **kwargs)
    return wrapper


class DiptoxPipeline:
    """Main processing class that coordinates various modules."""

    def __init__(self):
        self._check_initial_registration()
        self.chem_processor = ChemistryProcessor()
        self.data_handler = DataHandler()

        self.df: Optional[pd.DataFrame] = None
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

    def get_history(self) -> pd.DataFrame:
        """Returns the processing history as a DataFrame."""
        if not self._audit_log:
            return pd.DataFrame(columns=["Step", "Timestamp", "Rows Before", "Rows After", "Delta", "Details"])
        return pd.DataFrame(self._audit_log)

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
        df = self.data_handler.load_data(input_data=input_data, smiles_col=smiles_col, cas_col=cas_col,
                                         target_col=target_col, unit_col=unit_col, inchikey_col=inchikey_col,
                                         id_col=id_col, **kwargs)

        if 'Canonical SMILES' not in df.columns:
            df['Canonical SMILES'] = pd.Series(pd.NA, index=df.index, dtype="string")
        if 'Processing Log' not in df.columns:
            df['Processing Log'] = pd.Series(pd.NA, index=df.index, dtype="string")
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
            df['Is Valid'] = df['Is Valid'].astype("boolean")
        except Exception:
            df['Is Valid'] = pd.Series(False, index=df.index, dtype="boolean")

        self.df = df
        self.smiles_col = smiles_col
        self.cas_col = cas_col
        self.name_col = name_col
        self.target_col = target_col
        self.unit_col = unit_col
        self.inchikey_col = inchikey_col
        self.id_col = id_col
        self._units_standardized = False
        self._dedup_unit_settings = None

        source_name = input_data if isinstance(input_data, str) else "Memory/List"
        self._record_step("Data Loading", None, self.df, f"Source: {source_name}")

    def _ensure_dtypes_after_load(self) -> None:
        """
        Ensure critical columns use stable pandas extension dtypes to avoid
        FutureWarning when assigning incompatible values.
        """
        # string columns
        for col in ["Canonical SMILES", "Processing Log"]:
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
                   remove_mixtures: bool = False,
                   remove_inorganic: bool = True,
                   neutralize: bool = True,
                   reject_non_neutral: bool = False,
                   check_valid_atoms: bool = False,
                   strict_atom_check: bool = False,
                   remove_stereo: bool = False,
                   remove_isotopes: bool = True,
                   remove_hs: bool = True,
                   keep_largest_fragment: bool = True,
                   hac_threshold: int = 3,
                   sanitize: bool = True,
                   reject_radical_species: bool = True,
                   progress_callback: Optional[Callable[[int, int], None]] = None) -> pd.DataFrame:
        """
        Execute the chemical processing pipeline.\
        :param remove_salts: Whether to remove salts.
        :param remove_solvents: Whether to remove solvent molecules.
        :param remove_mixtures: Whether to remove mixtures.
        :param remove_inorganic: Whether to remove inorganic molecules.
        :param neutralize: Whether to neutralize charges.
        :param reject_non_neutral: Only retain the molecules whose formal charge is zero.
        :param check_valid_atoms: Whether to check for valid atoms.
        :param strict_atom_check: If True, remove the entire molecule if invalid atoms are found.
                                  If False, attempt to remove only the invalid atoms if they are not on the main chain.
        :param remove_stereo: Whether to remove stereochemistry.
        :param remove_isotopes: Whether to remove isotope information. Defaults to True.
        :param remove_hs: Whether to remove hydrogen atoms.
        :param keep_largest_fragment: Whether to keep the largest fragment.
        :param hac_threshold: Threshold for salt removal (heavy atoms count).
        :param sanitize: Whether to perform chemical sanitization.
        :param reject_radical_species: Molecules containing free radical atoms are directly rejected.
        :return: Processed DataFrame with results.
        """
        self._ensure_dtypes_after_load()
        df_start = self.df.copy()

        # Build processing steps
        steps: List[Callable[[Chem.Mol], Optional[Chem.Mol]]] = []
        step_descriptions: List[str] = []

        if remove_stereo:
            steps.append(self.chem_processor.remove_stereochemistry)
            step_descriptions.append("Stereo removal")

        if remove_isotopes:
            steps.append(self.chem_processor.remove_isotopes)
            step_descriptions.append("Isotope removal")

        if remove_hs:
            steps.append(self.chem_processor.remove_hydrogens)
            step_descriptions.append("Hydrogen removal")

        if remove_salts:
            salt_processor = partial(self.chem_processor.remove_salts)
            steps.append(salt_processor)
            step_descriptions.append("Salt removal")

        if remove_solvents:
            solvent_processor = partial(self.chem_processor.remove_solvents)
            steps.append(solvent_processor)
            step_descriptions.append("Solvent removal")

        if remove_mixtures:
            mixture_processor = partial(
                self.chem_processor.remove_mixtures,
                hac_threshold=hac_threshold,
                keep_largest=keep_largest_fragment
            )
            steps.append(mixture_processor)
            step_descriptions.append("Mixture removal")

        if remove_inorganic:
            steps.append(self.chem_processor.remove_inorganic)
            step_descriptions.append("Inorganic removal")

        if reject_radical_species:
            steps.append(self.chem_processor.reject_radicals)
            step_descriptions.append("Radical check")

        if neutralize:
            charge_processor = partial(
                self.chem_processor.neutralize_charges,
                reject_non_neutral=reject_non_neutral
            )
            steps.append(charge_processor)
            description = "Charge neutralization"
            if reject_non_neutral:
                description += " and Non-neutral rejection"
            step_descriptions.append(description)

        if check_valid_atoms:
            atom_validator = partial(self.chem_processor.effective_atom, strict=strict_atom_check)
            steps.append(atom_validator)
            step_descriptions.append("Atom validation")

        total_rows = len(self.df)

        # Process each molecule
        for idx, row in tqdm(self.df.iterrows(), total=len(self.df), desc="Processing"):
            smiles = row[self.smiles_col]
            comments = []

            mol = self.chem_processor.smiles_to_mol(smiles, sanitize)
            if mol is None:
                self._update_row(idx, False, "Invalid SMILES", None)
                continue

            # Execute processing steps
            for step, desc in zip(steps, step_descriptions):
                if mol is None:
                    break

                try:
                    processed = step(mol)
                    if processed is None:
                        comments.append(f"{desc} failed")
                        mol = None
                    else:
                        mol = processed
                except Exception as e:
                    comments.append(f"{desc} error: {str(e)}")
                    mol = None

            # Results handling
            if mol is not None:
                try:
                    std_smiles = self.chem_processor.standardize_smiles(mol)
                    self._update_row(idx, True, "; ".join(comments) if comments else "Success", std_smiles)
                except Exception as e:
                    self._update_row(idx, False, f"Standardization failed: {str(e)}", None)
            else:
                self._update_row(idx, False, "; ".join(comments), None)

            if progress_callback and idx % 10 == 0:
                progress_callback(idx + 1, total_rows)

        if progress_callback:
            progress_callback(total_rows, total_rows)
        self._preprocess_key = 1

        valid_count = self.df['Is Valid'].sum()
        invalid_count = total_rows - valid_count
        pipeline_order = " -> ".join(step_descriptions)
        self._record_step("Preprocessing", df_start, self.df,
                          f"Valid: {valid_count}, Invalid: {invalid_count}. Order: {pipeline_order}")
        return self.df

    def _update_row(self, idx, is_valid: bool, comment: str, smiles: Optional[str]) -> None:
        """Update the result row for a given index."""
        self.df.at[idx, 'Is Valid'] = bool(is_valid)
        self.df.at[idx, 'Processing Log'] = "" if comment is None else str(comment)
        self.df.at[idx, 'Canonical SMILES'] = (pd.NA if smiles is None else str(smiles))

    @check_data_loaded
    def standardize_units(self, standard_unit: Optional[str] = None,
                          conversion_rules: Optional[Dict[Tuple[str, str], str]] = None) -> None:
        """
        Orchestrates the standardization of units for the target column.
        :param standard_unit: The target unit to convert all values to.
        :param conversion_rules: A dictionary of conversion rules, e.g., {('mg/L', 'ug/L'): 'x * 1000'}.
        """
        df_start = self.df.copy()
        if not self.target_col or not self.unit_col:
            logger.info("Target column or unit column not specified, skipping unit standardization.")
            self._units_standardized = True
            return

        unique_units = [u for u in self.df[self.unit_col].dropna().unique() if u]
        if len(unique_units) <= 1:
            logger.info("Only one unit detected. No conversion necessary.")
            # Ensure a consistent '_new' column is created for the next step
            new_target_col = f"{self.target_col} (Standardized)"
            new_unit_col = f"{self.unit_col} (Standardized)"
            self.df[new_target_col] = self.df[self.target_col]
            self.df[new_unit_col] = self.df[self.unit_col]
            self.target_col = new_target_col
            self.unit_col = new_unit_col
            self._units_standardized = True
            return

        final_standard_unit = standard_unit
        is_gui_mode = os.environ.get("DIPTOX_GUI_MODE") == "true"

        if not final_standard_unit:
            if is_gui_mode:
                raise ValueError("A standard unit must be provided when multiple units exist.")
            final_standard_unit = self._select_standard_unit_interactively(unique_units)
            if not final_standard_unit:
                return

        unit_processor = UnitProcessor(rules=conversion_rules)

        rule_provider = None
        if not is_gui_mode:
            prompt_tracker = {'first_time': True}
            rule_provider = lambda from_unit, to_unit: self._get_rule_from_user(from_unit, to_unit, prompt_tracker)

        try:
            self.df, new_target_col, new_unit_col = unit_processor.standardize(
                df=self.df,
                target_col=self.target_col,
                unit_col=self.unit_col,
                standard_unit=final_standard_unit,
                rule_provider_callback=rule_provider
            )
            self.target_col = new_target_col
            self.unit_col = new_unit_col
            self._units_standardized = True
            self._record_step("Unit Standardization", df_start, self.df, f"Target: {final_standard_unit}")
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

    def config_deduplicator(self, condition_cols: Optional[List[str]] = None,
                            data_type: str = "continuous",
                            method: str = "auto",
                            priority: Optional[List[str]] = None,
                            p_threshold: float = 0.05,
                            custom_method: Optional[Callable] = None,
                            standard_unit: Optional[str] = None,
                            conversion_rules: Optional[Dict[Tuple[str, str], str]] = None,
                            log_transform: bool = False,
                            dropna_conditions: bool = False) -> None:
        """
        Configure the deduplicator device
        :param condition_cols: Data condition column (e.g. temperature, pressure, etc.)
        :param data_type: data type - discrete/continuous
        :param method: Existing method of data deduplication (e.g., auto, vote, priority, 3sigma, IQR.)
        :param priority: List of preferred values for discrete data
        :param p_threshold: Threshold of normal distribution
        :param custom_method: Custom method of data deduplication
        :param standard_unit: The target unit to standardize to before deduplication.
        :param conversion_rules: A dictionary of rules for unit conversion.
        :param log_transform: If True, applies a -log10 transformation to the target column.
        :param dropna_conditions: If True, drops rows with missing condition values. If False, groups them.
        """
        if standard_unit or conversion_rules:
            self._dedup_unit_settings = {'standard_unit': standard_unit, 'conversion_rules': conversion_rules}

        smiles_col = 'Canonical SMILES' if self._preprocess_key else self.smiles_col

        self.deduplicator = DataDeduplicator(
            smiles_col=smiles_col, target_col=self.target_col, condition_cols=condition_cols,
            data_type=data_type, method=method, p_threshold=p_threshold, priority=priority,
            custom_method=custom_method, log_transform=log_transform, dropna_conditions=dropna_conditions
        )

    @check_data_loaded
    def dataset_deduplicate(self, progress_callback: Optional[Callable] = None) -> None:
        """Execution deduplicator removal"""
        if not self.deduplicator:
            raise ValueError("Deduplicator not configured. Call config_deduplicator first.")
        df_start = self.df.copy()

        if self._dedup_unit_settings and not self._units_standardized:
            logger.info("Implicitly running unit standardization as part of deduplication.")
            self.standardize_units(
                standard_unit=self._dedup_unit_settings.get('standard_unit'),
                conversion_rules=self._dedup_unit_settings.get('conversion_rules')
            )
            df_start = self.df.copy()

        needs_standardization = False
        if self.target_col and self.unit_col and self.unit_col in self.df.columns:
            unique_count = self.df[self.unit_col].dropna().nunique()
            needs_standardization = unique_count > 1

        if needs_standardization and not self._units_standardized:
            raise ValueError(
                "Unit standardization is required but has not been performed. Please go to the 'Unit Standardization' step first.")

        self.deduplicator.target_col = self.target_col

        self.df = self.deduplicator.deduplicate(self.df, progress_callback=progress_callback)
        if self.target_col:
            self.target_col = self.target_col + "_new"

        method_name = self.deduplicator.method
        if self.deduplicator.log_transform:
            method_name += " (Log10 Transformed)"
        self._record_step("Deduplication", df_start, self.df, f"Method: {method_name}")

    @check_data_loaded
    def substructure_search(self, query_pattern: Union[str, List[str]],
                            is_smarts: bool = False) -> None:
        """
        Integrated search interface
        :param query_pattern: Molecular substructure
        :param is_smarts: Search mode (SMILES/SMARTS)
        """
        searcher = SubstructureSearcher(
            df=self.df,
            smiles_col='Canonical SMILES' if self._preprocess_key else self.smiles_col,
        )
        query_pattern_list = [query_pattern] if isinstance(query_pattern, str) else query_pattern
        for query_pattern in query_pattern_list:
            results = searcher.search(query_pattern, is_smarts)
            col_name = f'Substructure_{query_pattern}'
            if col_name not in self.df.columns:
                self.df[col_name] = pd.Series(False, index=self.df.index, dtype="boolean")
            else:
                try:
                    self.df[col_name] = self.df[col_name].astype("boolean")
                except Exception:
                    self.df[col_name] = pd.Series(False, index=self.df.index, dtype="boolean")

            for idx, _ in results['matches']:
                self.df.at[idx, col_name] = True

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

                pending_indices = [idx for idx in pending_indices if idx not in processed_indices]

            if pending_indices:
                self.df.loc[pending_indices, 'Query_Status'] = 'Failed'

            logger.info("Web request processing is complete.")

            if 'smiles' in request_set:
                if not self._preprocess_key:
                    self.smiles_col = 'smiles_from_web'
            if 'cas' in request:
                self.cas_col = 'cas_from_web'
            if 'name' in request:
                self.name_col = 'name_from_web'

            success_count = len(self.df[self.df['Query_Status'] == 'Success'])
            self._record_step("Web Request", None, self.df, f"Found data for {success_count} rows")

        except requests.exceptions.RequestException as e:
            logger.error(f"An error occurred on the network, and the web request was interrupted: {e}")
            self.df.loc[pending_indices, 'Query_Status'] = 'Network Error'
            return
        except Exception as e:
            logger.error(f"An unknown error occurred during the web request: {e}")
            self.df.loc[pending_indices, 'Query_Status'] = 'Unknown Error'
            raise

    @check_data_loaded
    def calculate_inchi(self) -> None:
        """
        Calculate InChI strings locally using RDKit based on the current SMILES column.
        No web request required.
        """
        smiles_col = 'Canonical SMILES' if self._preprocess_key else self.smiles_col
        inchi_col = 'InChI'

        if inchi_col not in self.df.columns:
            self.df[inchi_col] = pd.Series(pd.NA, index=self.df.index, dtype="string")

        logger.info("Calculating InChI from SMILES using RDKit...")

        count = 0
        for idx, row in tqdm(self.df.iterrows(), total=len(self.df), desc="InChI Calc"):
            smiles = row[smiles_col]
            if pd.isna(smiles):
                continue

            mol = self.chem_processor.smiles_to_mol(smiles, sanitize=True)
            inchi = self.chem_processor.mol_to_inchi(mol)

            if inchi:
                self.df.at[idx, inchi_col] = inchi
                count += 1

        logger.info(f"InChI calculation complete. Generated {count} InChI strings.")
        self._record_step("InChI Calculation", None, self.df, f"Calculated {count} InChIs")

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
        df_start = self.df.copy()
        if all(arg is None for arg in [min_heavy_atoms, max_heavy_atoms, min_total_atoms, max_total_atoms]):
            logger.warning("No filter criteria provided for filter_by_atom_count. No action taken.")
            return

        initial_count = len(self.df)
        smiles_col = 'Canonical SMILES' if self._preprocess_key else self.smiles_col

        def is_valid_by_atom_count(row):
            s = row[smiles_col]
            if pd.isna(s) or str(s).strip() == "":
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

        self.df = self.df[mask].reset_index(drop=True)
        final_count = len(self.df)
        logger.info(
            f"Filtered by atom count. Initial: {initial_count}, Final: {final_count}, Removed: {initial_count - final_count}")
        details = f"Heavy: {min_heavy_atoms}-{max_heavy_atoms}, Total: {min_total_atoms}-{max_total_atoms}"
        self._record_step("Filter Atom Count", df_start, self.df, details)

    @check_data_loaded
    def save_results(self, output_path: str, columns: Optional[List[str]] = None) -> None:
        """
        Save the processed results to a file.
        :param output_path: The output path where the results will be saved.
        :param columns: The columns to save (default saves all columns).
        """
        save_cols = columns if columns else self.df.columns.tolist()
        self.data_handler.save_data(self.df, output_path, save_cols, 'Canonical SMILES' if self._preprocess_key else self.smiles_col, self.id_col)

    # Chemical rule management interface
    def add_neutralization_rule(self, reactant: str, product: str) -> None:
        """
        Add a new neutralization rule to the list, ensuring the rule is valid and there are no conflicts.
        :param reactant: SMARTS string for the reactant.
        :param product: SMILES string for the product.
        """
        return self.chem_processor.add_neutralization_rule(reactant, product)

    def remove_neutralization_rule(self, reactant: str) -> None:
        """
        Remove a charge neutralization rule.
        :param reactant: SMARTS string for the reactant.
        """
        return self.chem_processor.remove_neutralization_rule(reactant)

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

    def display_processing_rules(self) -> None:
        """
        Displays the current chemical processing rules being used,
        including valid atoms, neutralization rules, salts, and solvents.
        """
        self.chem_processor.display_current_rules()
