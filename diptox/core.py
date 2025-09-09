# diptox/core.py
import pandas as pd
from typing import Optional, List, Union, Tuple, Callable, Dict
from functools import partial, wraps
from tqdm import tqdm
from rdkit import Chem
from .chem_processor import ChemistryProcessor
from .web_request import WebService
from .data_io import DataHandler
from .data_deduplicator import DataDeduplicator
from .substructure_search import SubstructureSearcher
from .logger import log_manager
logger = log_manager.get_logger(__name__)


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
        self.chem_processor = ChemistryProcessor()
        self.data_handler = DataHandler()

        self.df: Optional[pd.DataFrame] = None
        self.smiles_col: str = "Smiles"
        self.cas_col: Optional[str] = None
        self.target_col: Optional[str] = None
        self.inchikey_col = None
        self.id_col: Optional[str] = None

        self.deduplicator = None
        self.web_service = None
        self._preprocess_key = 0
        self.web_source = None

    def load_data(self,
                  input_data: Union[str, List[str], pd.DataFrame],
                  smiles_col: str = None,
                  cas_col: Optional[str] = None,
                  target_col: Optional[str] = None,
                  inchikey_col: Optional[str] = None,
                  id_col: Optional[str] = None,
                  **kwargs) -> None:
        """
        Load data and initialize columns for processing.
        :param input_data: Path to input data (.xlsx/.xls/.csv/.txt/.sdf/.smi), or a list, or a DataFrame.
        :param smiles_col: The column name containing SMILES strings.
        :param cas_col: The column name for CAS Numbers.
        :param target_col: The column name for target values (optional).
        :param inchikey_col: The column name for InChIKeys (optional).
        :param id_col: The column name for SMI file's SMILES ID (optional)
        :param sep: CSV file delimiter.
        """
        df = self.data_handler.load_data(input_data, smiles_col, cas_col, target_col, inchikey_col, id_col, **kwargs)
        processing_cols = {
            'Canonical smiles': None,
            'IsValid': False,
            'ProcessingComments': ''
        }
        cols_to_add = {k: v for k, v in processing_cols.items() if k not in df.columns}
        if cols_to_add:
            df = df.assign(**cols_to_add)
        self.df = df
        self.smiles_col = smiles_col
        self.cas_col = cas_col
        self.target_col = target_col
        self.inchikey_col = inchikey_col
        self.id_col = id_col

    @check_data_loaded
    def preprocess(self, remove_salts: bool = True,
                   remove_solvents: bool = True,
                   remove_mixtures: bool = False,
                   remove_inorganic: bool = True,
                   neutralize: bool = True,
                   check_valid_atoms: bool = False,
                   strict_atom_check: bool = False,
                   remove_stereo: bool = False,
                   remove_isotopes: bool = True,
                   remove_hs: bool = True,
                   keep_largest_fragment: bool = True,
                   hac_threshold: int = 3,
                   sanitize: bool = True) -> pd.DataFrame:
        """
        Execute the chemical processing pipeline.\
        :param remove_salts: Whether to remove salts.
        :param remove_solvents: Whether to remove solvent molecules.
        :param remove_mixtures: Whether to remove mixtures.
        :param remove_inorganic: Whether to remove inorganic molecules.
        :param neutralize: Whether to neutralize charges.
        :param check_valid_atoms: Whether to check for valid atoms.
        :param strict_atom_check: If True, remove the entire molecule if invalid atoms are found.
                                  If False, attempt to remove only the invalid atoms if they are not on the main chain.
        :param remove_stereo: Whether to remove stereochemistry.
        :param remove_isotopes: Whether to remove isotope information. Defaults to True.
        :param remove_hs: Whether to remove hydrogen atoms.
        :param keep_largest_fragment: Whether to keep the largest fragment.
        :param hac_threshold: Threshold for salt removal (heavy atoms count).
        :param sanitize: Whether to perform chemical sanitization.
        :return: Processed DataFrame with results.
        """
        # Build processing steps
        steps: List[Callable[[Chem.Mol], Optional[Chem.Mol]]] = []
        step_descriptions: List[str] = []

        if remove_stereo:
            steps.append(self.chem_processor.remove_stereochemistry)
            step_descriptions.append("Stereo removal")

        if remove_isotopes:
            steps.append(self.chem_processor.remove_isotopes)
            step_descriptions.append("Isotope removal")

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

        if neutralize:
            steps.append(self.chem_processor.neutralize_charges)
            step_descriptions.append("Charge neutralization")

        if check_valid_atoms:
            atom_validator = partial(self.chem_processor.effective_atom, strict=strict_atom_check)
            steps.append(atom_validator)
            step_descriptions.append("Atom validation")

        if remove_hs:
            steps.append(self.chem_processor.remove_hydrogens)
            step_descriptions.append("Hydrogen removal")

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

        self._preprocess_key = 1
        return self.df

    def _update_row(self, idx, is_valid: bool, comment: str, smiles: Optional[str]) -> None:
        """Update the result row for a given index."""
        self.df.at[idx, 'IsValid'] = is_valid
        self.df.at[idx, 'ProcessingComments'] = comment
        self.df.at[idx, 'Canonical smiles'] = smiles

    def config_deduplicator(self, condition_cols: Optional[List[str]] = None,
                            data_type: str = "continuous",
                            method: str = "auto",
                            p_threshold: float = 0.05,
                            custom_method: Optional[Callable] = None) -> None:
        """
        Configure the deduplicator device
        :param condition_cols: Data condition column (e.g. temperature, pressure, etc.)
        :param data_type: data type - discrete/continuous
        :param method: Existing method of data deduplication (e.g., auto, vote, 3sigma, IQR.)
        :param p_threshold: Threshold of normal distribution
        :param custom_method: Custom method of data deduplication
        """
        smiles_col = 'Canonical smiles' if self._preprocess_key else self.smiles_col

        self.deduplicator = DataDeduplicator(
            smiles_col=smiles_col,
            target_col=self.target_col,
            condition_cols=condition_cols,
            data_type=data_type,
            method=method,
            p_threshold=p_threshold,
            custom_method=custom_method
        )

    @check_data_loaded
    def data_deduplicate(self) -> None:
        """Execution deduplicator removal"""
        if not self.deduplicator:
            raise ValueError("Deduplicator not configured. Call config_deduplicator first.")

        self.df = self.deduplicator.deduplicate(self.df)
        if self.target_col:
            self.target_col = self.target_col + "_new"

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
            smiles_col='Canonical smiles' if self._preprocess_key else self.smiles_col,
        )
        query_pattern_list = [query_pattern] if isinstance(query_pattern, str) else query_pattern
        for query_pattern in query_pattern_list:
            results = searcher.search(query_pattern, is_smarts)
            self.df[f'Substructure_{query_pattern}'] = False
            for idx, _ in results['matches']:
                self.df.at[idx, f'Substructure_{query_pattern}'] = True

    def config_web_request(self, source: str = 'pubchem',
                           interval: int = 0.3,
                           retries: int = 3,
                           delay: int = 30,
                           max_workers: int = 4,
                           batch_limit: int = 1500,
                           rest_duration: int = 300,
                           chemspider_api_key: Optional[str] = None,
                           comptox_api_key: Optional[str] = None) -> None:
        """
        Initializes the WebService class
        :param source: Data source interface
        :param interval: Time interval in seconds
        :param retries: Number of retry attempts on failure.
        :param delay: Delay between retries (in seconds).
        :param max_workers: Maximum number of concurrent requests.
        :param batch_limit: Number of requests before taking a break.
        :param rest_duration: Duration of the break in seconds.
        :param chemspider_api_key: Chemspider API key.
        :param comptox_api_key: Comptox API key.
        """
        self.web_source = source
        self.web_service = WebService(source=source, interval=interval, retries=retries, delay=delay,
                                      max_workers=max_workers, batch_limit=batch_limit, rest_duration=rest_duration,
                                      chemspider_api_key=chemspider_api_key, comptox_api_key=comptox_api_key)

    @check_data_loaded
    def web_request(self, send: str, request: Union[str, List[str]]) -> None:
        """
        Add CAS numbers for valid molecules.
        :param send: What is used to request additional data? (smiles/cas)
        :param request: What identifier is requested? (smiles/cas/iupac)
        """
        send = send.strip().lower()
        if send not in ['smiles', 'cas', 'inchikey']:
            raise ValueError("Send must be 'smiles' or' cas' or 'inchikey'")
        request = [request] if isinstance(request, str) else request
        request = [prop.strip().lower() for prop in request]
        valid_props = {'smiles', 'cas', 'iupac', 'inchikey'}
        if invalid := set(request) - valid_props:
            raise ValueError(f"Invalid request properties: {invalid}")

        if send == 'cas':
            input_col = self.cas_col
            input_type = 'cas'
        elif send == 'smiles':
            input_col = 'Canonical smiles' if self._preprocess_key else self.smiles_col
            input_type = 'smiles'
        else:
            input_col = self.inchikey_col
            input_type = 'inchikey'
        if input_col not in self.df.columns:
            raise ValueError(f"Input column {input_col} not found in data.")

        results = self.web_service.get_properties(self.df[input_col].tolist(), request, input_type)
        column_map = {
            'cas': f"CAS Number_{self.web_source}",
            'iupac': f"IUPAC Name_{self.web_source}",
            'smiles': f"Smiles_{self.web_source}",
            'inchikey': f"InchiKey_{self.web_source}"
        }
        if 'smiles' in request:
            self.smiles_col, self._preprocess_key = column_map['smiles'], 0
        if 'cas' in request:
            self.cas_col = column_map['cas']
        if 'inchikey' in request:
            self.inchikey_col = column_map['inchikey']
        for prop in request:
            self.df[column_map[prop]] = [r.get(prop, None) for r in results]

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
        if all(arg is None for arg in [min_heavy_atoms, max_heavy_atoms, min_total_atoms, max_total_atoms]):
            logger.warning("No filter criteria provided for filter_by_atom_count. No action taken.")
            return

        initial_count = len(self.df)
        smiles_col = 'Canonical smiles' if self._preprocess_key else self.smiles_col

        def is_valid_by_atom_count(s):
            mol = self.chem_processor.smiles_to_mol(s, sanitize=False)  # Use internal method
            return self.chem_processor.validate_atom_count(
                mol,
                min_heavy_atoms,
                max_heavy_atoms,
                min_total_atoms,
                max_total_atoms
            )

        # Build the filter mask
        mask = self.df[smiles_col].apply(is_valid_by_atom_count)

        self.df = self.df[mask].reset_index(drop=True)
        final_count = len(self.df)
        logger.info(
            f"Filtered by atom count. Initial: {initial_count}, Final: {final_count}, Removed: {initial_count - final_count}")

    @check_data_loaded
    def save_results(self, output_path: str, columns: Optional[List[str]] = None) -> None:
        """
        Save the processed results to a file.
        :param output_path: The output path where the results will be saved.
        :param columns: The columns to save (default saves all columns).
        """
        save_cols = columns if columns else self.df.columns.tolist()
        self.data_handler.save_data(self.df, output_path, save_cols, 'Canonical smiles' if self._preprocess_key else self.smiles_col, self.id_col)

    # Chemical rule management interface
    def add_neutralization_rule(self, reactant: str, product: str) -> None:
        """
        Add a new neutralization rule to the list, ensuring the rule is valid and there are no conflicts.
        :param reactant: SMARTS string for the reactant.
        :param product: SMILES string for the product.
        """
        self.chem_processor.add_neutralization_rule(reactant, product)

    def remove_neutralization_rule(self, reactant: str) -> None:
        """
        Remove a charge neutralization rule.
        :param reactant: SMARTS string for the reactant.
        """
        self.chem_processor.remove_neutralization_rule(reactant)

    def manage_atom_rules(self, atoms: Union[str, List[str]], add: bool = True) -> None:
        """Manage atom validation rules."""
        atom_list = [atoms] if isinstance(atoms, str) else atoms
        for atom in atom_list:
            if add:
                self.chem_processor.add_effective_atom(atom)
            else:
                self.chem_processor.delete_effective_atom(atom)

    def manage_default_salt(self, salts: Union[str, List[str]], add: bool = True) -> None:
        """Manage salt validation rules."""
        salt_list = [salts] if isinstance(salts, str) else salts
        for salt in salt_list:
            if add:
                self.chem_processor.add_default_salt(salt)
            else:
                self.chem_processor.remove_default_salt(salt)

    def manage_default_solvent(self, solvents: Union[str, List[str]], add: bool = True) -> None:
        """Manage solvent validation rules."""
        solvent_list = [solvents] if isinstance(solvents, str) else solvents
        for solvent in solvent_list:
            if add:
                self.chem_processor.add_default_solvents(solvent)
            else:
                self.chem_processor.remove_default_solvents(solvent)

    def display_processing_rules(self) -> None:
        """
        Displays the current chemical processing rules being used,
        including valid atoms, neutralization rules, salts, and solvents.
        """
        self.chem_processor.display_current_rules()
