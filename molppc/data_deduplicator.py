# molppc/data_deduplicator.py
import pandas as pd
import numpy as np
from typing import List, Optional, Dict, Any, Callable, Tuple
from .logger import log_manager
logger = log_manager.get_logger(__name__)


class DataDeduplicator:
    """Data Deduplication Processor"""

    def __init__(self, smiles_col: str = "Smiles",
                 target_col: Optional[str] = None,
                 condition_cols: Optional[List[str]] = None,
                 data_type: str = "continuous",
                 method: str = "auto",
                 custom_method: Optional[Callable[[pd.Series], Tuple[pd.Series, str]]] = None):
        """
        :param smiles_col: The name of the SMILES column
        :param target_col: The name of the target value column (optional)
        :param condition_cols: Columns representing conditions (e.g., temperature, pressure, etc.)
        :param data_type: Data type - "discrete" or "continuous"
        :param method: Existing method of data deduplication (e.g., auto, vote, 3sigma, IQR.)
        :param custom_method: Custom method of data deduplication
        """
        self.smiles_col = smiles_col
        self.target_col = target_col
        self.condition_cols = condition_cols or []
        self.data_type = data_type
        self.method = method
        self.custom_method = custom_method

        if custom_method and not callable(custom_method):
            raise ValueError("custom_outlier_handler must be a callable function")

        if data_type not in ["discrete", "continuous", None]:
            raise ValueError("Invalid data_type. Must be 'discrete' or 'continuous'")

        if target_col and not condition_cols:
            logger.warning("Target column provided but no condition columns specified")

    def deduplicate(self, df: pd.DataFrame) -> pd.DataFrame:
        """Main deduplication method"""
        self._validate_columns(df)
        df = df[df[self.smiles_col].notna()]

        group_keys = [self.smiles_col] + self.condition_cols

        grouped = df.groupby(group_keys, group_keys=False, sort=False)

        if self.target_col:
            return self._process_with_target(grouped)
        return self._process_without_target(grouped)

    def _validate_columns(self, df: pd.DataFrame):
        """Validate column existence"""
        required_cols = [self.smiles_col]
        if self.target_col:
            required_cols.append(self.target_col)

        missing = [col for col in required_cols if col not in df.columns]
        if missing:
            raise KeyError(f"Missing required columns: {missing}")

    def _process_without_target(self, grouped) -> pd.DataFrame:
        """Simple deduplication without target value"""
        logger.info("Performing simple deduplication by SMILES")
        return grouped.first().reset_index()

    def _process_with_target(self, grouped) -> pd.DataFrame:
        """Complex deduplication with target value"""
        logger.info(f"Processing deduplication with target ({self.data_type} data)")

        processed = []
        for name, group in grouped:
            if len(group) == 1:
                processed.append(self._mark_record(group, method='no change'))
                continue

            if self.data_type == "discrete":
                processed_group = self._handle_discrete(group)
            else:
                processed_group = self._handle_continuous(group)

            if not processed_group.empty:
                processed.append(processed_group)

        return pd.concat(processed).reset_index(drop=True)

    def _handle_discrete(self, group: pd.DataFrame) -> pd.DataFrame:
        """Handle discrete data"""
        counts = group[self.target_col].value_counts()
        max_count = counts.max()

        if len(counts[counts == max_count]) > 1:
            logger.debug(f"Tie detected in group: {group[self.smiles_col].iloc[0]}")
            return pd.DataFrame()

        selected = counts.idxmax()
        return self._mark_record(group[group[self.target_col] == selected].head(1), method='vote')

    def _handle_continuous(self, group: pd.DataFrame) -> pd.DataFrame:
        """Handle continuous data"""
        n = len(group)
        values = group[self.target_col]

        if n <= 3:
            return self._handle_small_group(values, group)

        if self.custom_method:
            clean_values, method = self.custom_method(values)
        else:
            clean_values, method = self._remove_outliers(values)

        final_value = clean_values.mean()
        final_record = group.iloc[[values.sub(final_value).abs().argmin()]].copy()
        final_record[self.target_col] = final_value
        return self._mark_record(final_record, method=method)

    def _handle_small_group(self, values: pd.Series, group: pd.DataFrame) -> pd.DataFrame:
        """Handle small sample groups"""
        if len(values) == 1:
            return self._mark_record(group, method="no change")

        final_value = values.mean()
        final_record = group.iloc[[values.sub(final_value).abs().argmin()]].copy()
        final_record[self.target_col] = final_value
        return self._mark_record(final_record, method="<=3(mean)")

    def _remove_outliers(self, values: pd.Series):
        """Outlier removal (3sigma/IQR)"""
        if self.data_type == "continuous":
            if self.method == "auto":
                # Automatically select method: use IQR for non-normal distributions
                if self._is_normal_distribution(values):
                    return self._3sigma_filter(values), '3sigma'
                return self._iqr_filter(values), 'IQR'
            elif self.method == "IQR":
                return self._iqr_filter(values), 'IQR'
            elif self.method == "3sigma":
                return self._3sigma_filter(values), '3sigma'

    @staticmethod
    def _is_normal_distribution(values: pd.Series, p_threshold: float = 0.05) -> bool:
        """Normal distribution test (Shapiro-Wilk)"""
        from scipy.stats import shapiro
        stat, p = shapiro(values)
        return p > p_threshold

    @staticmethod
    def _3sigma_filter(values: pd.Series, n_sigma: int = 3) -> pd.Series:
        """3-sigma filtering"""
        mean = values.mean()
        std = values.std()
        lower = mean - n_sigma * std
        upper = mean + n_sigma * std
        return values[(values >= lower) & (values <= upper)]

    @staticmethod
    def _iqr_filter(values: pd.Series, k: float = 1.5) -> pd.Series:
        """Interquartile range (IQR) filtering"""
        q1 = values.quantile(0.25)
        q3 = values.quantile(0.75)
        iqr = q3 - q1
        lower = q1 - k * iqr
        upper = q3 + k * iqr
        return values[(values >= lower) & (values <= upper)]

    def _mark_record(self, record: pd.DataFrame,
                     method: Optional[str] = None) -> pd.DataFrame:
        """Mark processed record"""
        record = record.copy()
        if method:
            record["DeduplicateMethod"] = method
        return record

    @classmethod
    def create_pipeline(cls, steps: List[Dict[str, Any]]) -> List['DataDeduplicator']:
        """Create a processing pipeline"""
        return [cls(**step_config) for step_config in steps]
