# diptox/data_deduplicator.py
import pandas as pd
import numpy as np
from typing import List, Optional, Dict, Any, Callable, Tuple
import warnings
from .logger import log_manager
logger = log_manager.get_logger(__name__)


class DataDeduplicator:
    """Data Deduplication Processor"""

    def __init__(self, smiles_col: str = "Smiles",
                 target_col: Optional[str] = None,
                 condition_cols: Optional[List[str]] = None,
                 data_type: str = "continuous",
                 method: str = "auto",
                 p_threshold: float = 0.05,
                 priority: Optional[List[str]] = None,
                 custom_method: Optional[Callable[[pd.Series], Tuple[pd.Series, str]]] = None,
                 log_transform: bool = False,
                 dropna_conditions: bool = False):
        """
        :param smiles_col: The name of the SMILES column
        :param target_col: The name of the target value column (optional)
        :param condition_cols: Columns representing conditions (e.g., temperature, pressure, etc.)
        :param data_type: Data type - "discrete" or "continuous"
        :param method: Existing method of data deduplication (e.g., auto, vote, 3sigma, IQR.)
        :param p_threshold: Threshold of normal distribution
        :param priority: List of values in descending order of priority for discrete deduplication
        :param custom_method: Custom method of data deduplication
        :param log_transform: If True, applies a -log10 transformation to the target column before continuous deduplication.
        :param dropna_conditions: If True, rows with NaN in condition columns are dropped.
        """
        self.smiles_col = smiles_col
        self.target_col = target_col
        self.condition_cols = condition_cols or []
        self.data_type = data_type
        self.method = method
        self.custom_method = custom_method
        self._p_threshold = p_threshold
        self.priority_list = priority
        self.log_transform = log_transform
        self.dropna_conditions = dropna_conditions

        if custom_method and not callable(custom_method):
            raise ValueError("custom_outlier_handler must be a callable function")

        if data_type not in ["smiles", "discrete", "continuous", None]:
            raise ValueError("Invalid data_type. Must be 'discrete', 'continuous', and 'smiles'")

        if target_col and not condition_cols:
            logger.warning("Target column provided but no condition columns specified")

        if self.target_col and self.data_type == 'smiles':
            logger.info(
                f"Data type is 'smiles'. The target column '{self.target_col}' will be ignored during deduplication.")

    def deduplicate(self, df: pd.DataFrame, progress_callback: Optional[Callable] = None) -> pd.DataFrame:
        """Main deduplication method"""
        self._validate_columns(df)
        df = df[df[self.smiles_col].notna()].copy()

        if self.target_col and self.data_type == 'continuous' and self.log_transform:
            logger.info(f"Applying -log10 transformation to the target column '{self.target_col}'.")
            initial_rows = len(df)
            positive_mask = pd.to_numeric(df[self.target_col], errors='coerce') > 0
            if not positive_mask.all():
                df = df[positive_mask]
                removed_count = initial_rows - len(df)
                logger.warning(
                    f"Removed {removed_count} rows with non-positive target values before log transformation.")
            df[self.target_col] = -np.log10(pd.to_numeric(df[self.target_col], errors='coerce'))

        group_keys = [self.smiles_col] + self.condition_cols
        grouped = df.groupby(group_keys, group_keys=False, sort=False, dropna=self.dropna_conditions)

        if self.data_type == 'smiles':
            return self._process_without_target(grouped)

        if self.target_col:
            return self._process_with_target(grouped, self._p_threshold, progress_callback)
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
        deduplicated_df = grouped.first().reset_index()
        return deduplicated_df.assign(**{'Deduplication Strategy': 'smiles_only'})

    def _process_with_target(self, grouped, p_threshold: float,
                             progress_callback: Optional[Callable] = None) -> pd.DataFrame:
        """Complex deduplication with target value"""
        logger.info(f"Processing deduplication with target ({self.data_type} data)")

        processed = []
        total_groups = len(grouped)
        for i, (name, group) in enumerate(grouped):
            if progress_callback and i % 50 == 0:
                progress_callback(i + 1, total_groups)
            valid_group = group.dropna(subset=[self.target_col])

            if valid_group.empty:
                logger.debug(f"Dropped group {name} due to all NaN targets.")
                continue

            if len(valid_group) == 1:
                valid_group = valid_group.copy()
                valid_group[self.target_col + '_new'] = valid_group[self.target_col]
                processed.append(self._mark_record(valid_group, method='no change'))
                continue

            if self.data_type == "discrete":
                processed_group = self._handle_discrete(valid_group)
            else:
                processed_group = self._handle_continuous(valid_group, p_threshold)

            if not processed_group.empty:
                processed.append(processed_group)

        if progress_callback:
            progress_callback(total_groups, total_groups)

        if not processed:
            original_columns = list(grouped.obj.columns)
            new_columns = [self.target_col + '_new', 'Deduplication Strategy']
            final_columns = original_columns + new_columns
            return pd.DataFrame(columns=final_columns)

        return pd.concat(processed).reset_index(drop=True)

    def _handle_discrete(self, group: pd.DataFrame) -> pd.DataFrame:
        """Handle discrete data"""
        if self.priority_list:
            for priority_val in self.priority_list:
                matches = group[group[self.target_col] == priority_val]
                if matches.empty:
                    matches = group[group[self.target_col].astype(str) == str(priority_val)]
                if not matches.empty:
                    final_record = matches.head(1).copy()
                    final_record[self.target_col + '_new'] = final_record[self.target_col]
                    return self._mark_record(final_record, method=f'priority_{priority_val}')

        counts = group[self.target_col].value_counts()
        max_count = counts.max()

        if len(counts[counts == max_count]) > 1:
            logger.debug(f"Tie detected in group: {group[self.smiles_col].iloc[0]}")
            return pd.DataFrame()

        selected = counts.idxmax()
        final_record = group[group[self.target_col] == selected].head(1).copy()
        final_record[self.target_col + '_new'] = selected
        return self._mark_record(final_record, method='vote')

    def _handle_continuous(self, group: pd.DataFrame, p_threshold: float) -> pd.DataFrame:
        """Handle continuous data"""
        n = len(group)
        values = group[self.target_col]

        if n <= 3:
            return self._handle_small_group(values, group)

        if self.custom_method:
            clean_values, method = self.custom_method(values)
        else:
            clean_values, method = self._remove_outliers(values, p_threshold)

        if clean_values.empty:
            logger.warning(
                f"All values in group for SMILES '{group[self.smiles_col].iloc[0]}' "
                f"were removed as outliers. Skipping this group."
            )
            return pd.DataFrame()

        final_value = clean_values.mean()
        valid_indices = clean_values.index
        valid_group = group.loc[valid_indices]
        best_idx = valid_group[self.target_col].sub(final_value).abs().idxmin()
        final_record = group.loc[[best_idx]].copy()
        final_record[self.target_col + '_new'] = final_value
        return self._mark_record(final_record, method=method)

    def _handle_small_group(self, values: pd.Series, group: pd.DataFrame) -> pd.DataFrame:
        """Handle small sample groups"""
        # if len(values) == 1:
        #     group[self.target_col + '_new'] = group[self.target_col]
        #     return self._mark_record(group, method="no change")

        final_value = values.mean()
        final_record = group.iloc[[values.sub(final_value).abs().argmin()]].copy()
        final_record[self.target_col + '_new'] = final_value
        return self._mark_record(final_record, method="<=3(mean)")

    def _remove_outliers(self, values: pd.Series, p_threshold: float):
        """Outlier removal (3sigma/IQR)"""
        if self.data_type == "continuous":
            if self.method == "auto":
                # Automatically select method: use IQR for non-normal distributions
                if self._is_normal_distribution(values, p_threshold):
                    return self._3sigma_filter(values), '3sigma'
                return self._iqr_filter(values), 'IQR'
            elif self.method == "IQR":
                return self._iqr_filter(values), 'IQR'
            elif self.method == "3sigma":
                return self._3sigma_filter(values), '3sigma'

    @staticmethod
    def _is_normal_distribution(values: pd.Series, p_threshold: float) -> bool:
        """Normal distribution test (Shapiro-Wilk)"""
        from scipy.stats import shapiro
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Input data for shapiro has range zero.*", category=UserWarning)
            try:
                stat, p = shapiro(values)
            except ValueError:
                return True
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
            record["Deduplication Strategy"] = method
        return record

    @classmethod
    def create_pipeline(cls, steps: List[Dict[str, Any]]) -> List['DataDeduplicator']:
        """Create a processing pipeline"""
        return [cls(**step_config) for step_config in steps]
