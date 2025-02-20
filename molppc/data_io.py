# molppc/data_io.py
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
                  target_col: Optional[str] = None) -> pd.DataFrame:
        """
        Unified data loading entry point
        :param input_data: Supports three i nput types:
            - File path (csv/xlsx/xls/txt)
            - SMILES list
            - Pre-loaded DataFrame
        :param smiles_col: SMILES column name
        :param cas_col: The column name for CAS Numbers.
        :param target_col: Target value column name (optional)
        """
        if isinstance(input_data, str):
            return DataHandler._load_from_file(input_data, smiles_col, cas_col, target_col)
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
                        target_col: Optional[str] = None) -> pd.DataFrame:
        """Load data from file"""
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File does not exist: {file_path}")

        if file_path.endswith('.csv'):
            df = pd.read_csv(file_path)
        elif file_path.endswith(('.xls', '.xlsx')):
            df = pd.read_excel(file_path)
        elif file_path.endswith('.txt'):
            df = pd.read_csv(file_path, sep='\t')
        else:
            logger.error("Only the .csv/.xlsx/.xls/.txt file format is supported")

        if smiles_col is not None and smiles_col not in df.columns:
            raise logger.error(f"SMILES column '{smiles_col}' does not exist in the file")
        if cas_col is not None and cas_col not in df.columns:
            raise logger.error(f"CAS column '{cas_col}' does not exist in the file")
        if target_col is not None and target_col not in df.columns:
            raise logger.error(f"The target value column '{target_col}' does not exist in the file")

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
    def save_data(df: pd.DataFrame, output_path: str, columns: list):
        """
        Save processing results
        :param df: The DataFrame that contains the result
        :param output_path: Output path of the data
        :param columns: Columns to save
        """
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
                    df[columns].to_csv(output_path, index=False, sep='\t')
                else:
                    logger.warning(f"Unsupported file format. The file will be saved as csv by default.")
                    output_path += '.csv'
                    df[columns].to_csv(output_path, index=False, encoding='utf-8')
                logger.info(f"File saved successfully: {current_output_path}")
                break
            except (PermissionError, IOError, OSError) as e:
                logger.error(f"Unable to save file: {str(e)}")
                choice = input("Do you want to save again? (Y/N): ").strip().lower()
                if choice == 'y':
                    continue
                else:
                    logger.warning("The user canceled the save operation")
                    break
            except Exception as e:
                logger.error(f"Unknown error: {str(e)}")
                break
