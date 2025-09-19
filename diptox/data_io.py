# diptox/data_io.py
import pandas as pd
from typing import Union, List, Optional, Dict, Any
import os
from .logger import log_manager
logger = log_manager.get_logger(__name__)


class DataHandler:
    """Data loading and saving"""

    @staticmethod
    def load_data(input_data: Union[str, List[str], pd.DataFrame],
                  smiles_col: str = None,
                  cas_col: Optional[str] = None,
                  name_col: Optional[str] = None,
                  target_col: Optional[str] = None,
                  inchikey_col: Optional[str] = None,
                  id_col: Optional[str] = None,
                  **kwargs) -> pd.DataFrame:
        """
        Unified data loading entry point
        :param input_data: Supports three i nput types:
            - File path (csv/xlsx/xls/txt)
            - SMILES list
            - Pre-loaded DataFrame
        :param smiles_col: SMILES column name
        :param cas_col: The column name for CAS Numbers.
        :param name_col: Name column name
        :param target_col: Target value column name (optional)
        :param inchikey_col: Inchikey column name (optional)
        :param id_col: SMI file's SMILES ID column name (optional)
        """
        if isinstance(input_data, str):
            return DataHandler._load_from_file(input_data, smiles_col, cas_col, name_col, target_col, inchikey_col, id_col, **kwargs)
        elif isinstance(input_data, list):
            return DataHandler._load_from_list(input_data, smiles_col)
        elif isinstance(input_data, pd.DataFrame):
            return DataHandler._load_from_dataframe(input_data, smiles_col, target_col)
        else:
            logger.error(f"Unsupported input types: {type(input_data)}")

    @staticmethod
    def _load_from_file(file_path: str,
                        smiles_col: Optional[str] = None,
                        cas_col: Optional[str] = None,
                        name_col: Optional[str] = None,
                        target_col: Optional[str] = None,
                        inchikey_col: Optional[str] = None,
                        id_col: Optional[str] = None,
                        **kwargs) -> pd.DataFrame:
        """Load data from file"""
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File does not exist: {file_path}")

        if file_path.endswith('.csv'):
            df = pd.read_csv(file_path, **kwargs)
        elif file_path.endswith(('.xls', '.xlsx')):
            if 'sheet_name' in kwargs:
                df = pd.read_excel(file_path, **kwargs)
            else:
                xls = pd.ExcelFile(file_path)
                sheet_names = xls.sheet_names
                if len(sheet_names) == 1:
                    df = pd.read_excel(file_path, sheet_name=sheet_names[0], **kwargs)
                else:
                    print("The file contains multiple sheets. Please select one:")
                    for i, sheet_name in enumerate(sheet_names):
                        print(f"{i + 1}: {sheet_name}")

                    while True:
                        user_input = input("Enter the sheet number or name: ").strip()
                        try:
                            choice = int(user_input)
                            if 1 <= choice <= len(sheet_names):
                                df = pd.read_excel(file_path, sheet_name=sheet_names[choice - 1], **kwargs)
                                break
                            else:
                                print(f"Invalid number. Please enter a number between 1 and {len(sheet_names)}.")
                        except ValueError:
                            if user_input in sheet_names:
                                df = pd.read_excel(file_path, sheet_name=user_input)
                                break
                            else:
                                print(f"Invalid sheet name. Please enter a valid sheet name or number.")
        elif file_path.endswith('.txt'):
            df = pd.read_csv(file_path, sep='\t', **kwargs)
        elif file_path.endswith('.sdf'):
            df = DataHandler._load_sdf(file_path=file_path, smiles_col=smiles_col, target_col=target_col)
        elif file_path.endswith('.smi'):
            df = DataHandler._load_smi(file_path=file_path, smiles_col=smiles_col, id_col=id_col, **kwargs)
        else:
            logger.error("Only the .csv/.xlsx/.xls/.txt file format is supported")

        if smiles_col is not None and smiles_col not in df.columns:
            logger.error(f"SMILES column '{smiles_col}' does not exist in the file")
            raise
        if cas_col is not None and cas_col not in df.columns:
            logger.error(f"CAS column '{cas_col}' does not exist in the file")
            raise
        if inchikey_col is not None and inchikey_col not in df.columns:
            logger.error(f"Inchikey column '{inchikey_col}' does not exist in the file")
            raise
        if target_col is not None and target_col not in df.columns:
            logger.error(f"The target value column '{target_col}' does not exist in the file")
            raise
        if name_col is not None and name_col not in df.columns:
            logger.error(f"The name column '{name_col}' does not exist in the file")
            raise

        return df

    @staticmethod
    def _load_from_list(smiles_list: List[str],
                        smiles_col: str) -> pd.DataFrame:
        """Load data from a list directly"""
        if not all(isinstance(s, str) for s in smiles_list):
            logger.error("The SMILES list must all be of string type")

        if smiles_col is None:
            logger.warning("The SMILES column name are not recommended to be left blank")
        data = {smiles_col: smiles_list}

        return pd.DataFrame(data)

    @staticmethod
    def _load_from_dataframe(df: pd.DataFrame,
                             smiles_col: str,
                             target_col: Optional[str] = None) -> pd.DataFrame:
        """Load data from an existing DataFrame"""
        required_cols = [smiles_col]
        if target_col:
            required_cols.append(target_col)

        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            logger.error(f"Missing necessary columns: {missing_cols}")

        return df.copy()

    @staticmethod
    def _load_sdf(file_path: str, smiles_col: str, target_col: str) -> pd.DataFrame:
        try:
            from rdkit import Chem
            from rdkit.Chem import PandasTools
            df = PandasTools.LoadSDF(file_path, removeHs=False, sanitize=False)
            original_columns = df.columns.tolist()
            if smiles_col:
                if smiles_col not in df.columns:
                    available = [c for c in original_columns if 'smiles' in c]
                    logger.error(f"'{smiles_col}' field not found in SDF, SMILES field available: {available}")
                    raise

            if target_col and target_col not in df.columns:
                logger.warning(f"target column '{target_col}' does not exist, available fields: {original_columns}")

            return df.rename(columns={smiles_col.lower(): 'smiles'})
        except Exception as e:
            logger.error(f"The SDF file fails to be parsed: {str(e)}")

    @staticmethod
    def _load_smi(file_path: str,
                  smiles_col: str = 'smiles',
                  id_col: Optional[str] = None,
                  **kwargs) -> pd.DataFrame:
        try:
            with open(file_path) as f:
                first_line = f.readline()
            sep = '\t' if '\t' in first_line else ' '

            if smiles_col in first_line:
                header = 0
            else:
                header = None

            df = pd.read_csv(file_path,
                             sep=sep,
                             header=header,
                             names=[col for col in [smiles_col, id_col] if col is not None][:len(first_line.split(sep))],
                             **kwargs)

            df[smiles_col] = df[smiles_col].str.strip()
            return df

        except Exception as e:
            logger.error(f"Parsing the SMI file failed: {str(e)}")
            raise

    @staticmethod
    def save_data(df: pd.DataFrame, output_path: str, columns: list, smiles_col: str, id_col: Optional[str] = None,):
        """Save processing results."""
        missing_cols = [col for col in columns if col not in df.columns]
        if missing_cols:
            logger.error(f"Columns {missing_cols} not found in DataFrame")

        directory = os.path.dirname(output_path)
        if directory:
            os.makedirs(directory, exist_ok=True)

        original_output_path = output_path
        while True:
            try:
                current_output_path = original_output_path
                if output_path.endswith('.csv'):
                    df[columns].to_csv(output_path, index=False, encoding='utf-8')
                elif output_path.endswith(('.xls', '.xlsx')):
                    df[columns].to_excel(output_path, index=False)
                elif output_path.endswith('.txt'):
                    df[columns].to_csv(output_path, index=False, sep='\t', encoding='utf-8')
                elif output_path.endswith('.sdf'):
                    from rdkit import Chem
                    from rdkit.Chem import PandasTools
                    if 'ROMol' not in df.columns:
                        df['ROMol'] = df[smiles_col].apply(
                            lambda s: Chem.MolFromSmiles(str(s)))
                    other = [c for c in columns if c != 'ROMol']
                    columns = ['ROMol'] + other
                    PandasTools.WriteSDF(df[columns], output_path, molColName='ROMol', properties=other)
                elif output_path.endswith('.smi'):
                    if id_col is None:
                        df[smiles_col].to_csv(output_path, sep='\t', header=True, index=False)
                    else:
                        df[[smiles_col, id_col]].to_csv(output_path, sep='\t', header=True, index=False)
                else:
                    logger.warning(f"Unsupported file format. The file will be saved as csv by default.")
                    output_path += '.csv'
                    df[columns].to_csv(output_path, index=False, encoding='utf-8')
                logger.info(f"File saved successfully: {current_output_path}")
                break
            except (PermissionError, IOError, OSError) as e:
                logger.error(f"Unable to save file: {str(e)}")
                choice = input("Do you want to save again? (Y/N): ").strip().lower()
                if choice in {'y', 'yes'}:
                    continue
                else:
                    logger.warning("The user canceled the save operation")
                    break
            except Exception as e:
                logger.error(f"Unknown error: {str(e)}")
                break
