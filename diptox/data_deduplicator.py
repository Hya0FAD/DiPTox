# diptox/data_deduplicator.py
import pandas as pd
import numpy as np
from typing import List, Optional, Dict, Any, Callable, Tuple
import warnings
from copy import deepcopy
from .logger import log_manager
from .unit_processor import TARGET_SCALES_ATTR, get_target_scale
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
                 dropna_conditions: bool = False,
                 aggregation: str = "mean"):
        """
        :param smiles_col: The name of the SMILES column
        :param target_col: The name of the target value column (optional)
        :param condition_cols: Columns representing conditions (e.g., temperature, pressure, etc.)
        :param data_type: Data type - "discrete" or "continuous"
        :param method: Continuous outlier filtering: auto, 3sigma or IQR.
                       Discrete selection: vote or priority.
        :param aggregation: Continuous aggregation after filtering: mean, max or min.
                            Groups of <=3 skip built-in filtering. Extrema ties keep the first row.
        :param p_threshold: Threshold of normal distribution
        :param priority: Required ordered values for discrete method='priority'; ignored by vote
                         and other methods. If no priority value occurs, fall back to voting.
        :param custom_method: Custom method of data deduplication
        :param log_transform: If True, applies a -log10 transformation to the target column before continuous deduplication.
        :param dropna_conditions: If True, rows with NaN in condition columns are dropped.
        """
        self.smiles_col = smiles_col
        self.target_col = target_col
        self.condition_cols = condition_cols or []
        self.data_type = data_type
        self.method = method
        self.aggregation = aggregation
        self.custom_method = custom_method
        self._p_threshold = p_threshold
        self.priority_list = priority if data_type == 'discrete' and method == 'priority' else None
        self.log_transform = log_transform
        self.dropna_conditions = dropna_conditions
        self.exclusion_reasons: Dict[Any, str] = {}
        # Structured membership for callers that maintain private row identities.
        # It follows the public source-row audit: all group inputs, including outliers.
        self.source_indices: Dict[Any, Tuple[Any, ...]] = {}

        if custom_method and not callable(custom_method):
            raise ValueError("custom_outlier_handler must be a callable function")

        if data_type not in ["smiles", "discrete", "continuous", None]:
            raise ValueError("Invalid data_type. Must be 'discrete', 'continuous', and 'smiles'")

        if aggregation not in {'mean', 'max', 'min'}:
            raise ValueError("aggregation must be 'mean', 'max' or 'min'")
        if data_type == 'continuous' and method not in {'auto', '3sigma', 'IQR'}:
            raise ValueError("Continuous method must be 'auto', '3sigma' or 'IQR'; use aggregation for 'mean', 'max' or 'min'")
        if data_type not in {'continuous', None} and aggregation != 'mean':
            raise ValueError("aggregation requires continuous data")

        if method == 'priority' and (data_type != 'discrete' or not self.priority_list):
            raise ValueError("method='priority' requires discrete data and a non-empty priority list")

        if target_col and not condition_cols:
            logger.info("No condition columns specified; standardized structure is the sole grouping key.")

        if self.target_col and self.data_type == 'smiles':
            logger.info(
                f"Data type is 'smiles'. The target column '{self.target_col}' will be ignored during deduplication.")

    def deduplicate(self, df: pd.DataFrame, progress_callback: Optional[Callable] = None) -> pd.DataFrame:
        """Main deduplication method"""
        self.exclusion_reasons = {}
        self.source_indices = {}
        self._validate_columns(df)
        if not df.index.is_unique:
            raise ValueError(
                "DataDeduplicator requires a unique DataFrame index. Use "
                "DiptoxPipeline.load_data to assign private unique row IDs automatically."
            )
        source_df = df
        df = df.copy()

        missing_structure = df[self.smiles_col].isna() | df[self.smiles_col].astype("string").str.strip().eq("")
        self._mark_excluded(df.index[missing_structure], "Missing standardized structure")
        df = df[~missing_structure].copy()

        if self.target_col and self.data_type == 'continuous':
            numeric_target = pd.to_numeric(df[self.target_col], errors='coerce')
            nonfinite_target = numeric_target.notna() & ~np.isfinite(numeric_target)
            invalid_target = numeric_target.isna() | nonfinite_target
            self._mark_excluded(df.index[nonfinite_target], "Non-finite continuous target")
            self._mark_excluded(
                df.index[invalid_target],
                "Invalid or missing continuous target",
            )
            dropped = int(invalid_target.sum())
            df = df[~invalid_target].copy()
            df[self.target_col] = numeric_target.loc[df.index]
            if dropped > 0:
                logger.warning(f"From the column '{self.target_col}', {dropped} records containing invalid, missing, or non-finite values have been removed.")
        elif self.target_col and self.data_type != 'smiles':
            missing_target = df[self.target_col].isna()
            self._mark_excluded(df.index[missing_target], "Missing target value")
            df = df[~missing_target].copy()

        transform_mode = self.log_transform
        if isinstance(transform_mode, bool):
            transform_mode = "-log10" if transform_mode else "None"
        if transform_mode is None:
            transform_mode = "None"
        if transform_mode not in {'None', 'log10', '-log10'}:
            raise ValueError("log_transform must be 'None', 'log10', '-log10', or a boolean")

        output_scale = get_target_scale(source_df, self.target_col)
        apply_transform = self.target_col and self.data_type == 'continuous' and transform_mode != 'None'
        if apply_transform and output_scale != 'linear':
            if output_scale != transform_mode:
                raise ValueError(
                    f"Target '{self.target_col}' is already on the {output_scale} scale; "
                    f"cannot apply {transform_mode}. Use the original linear target."
                )
            apply_transform = False

        if apply_transform:
            logger.info(f"Applying {transform_mode} transformation to the target column '{self.target_col}'.")
            initial_rows = len(df)
            positive_mask = df[self.target_col] > 0

            if not positive_mask.all():
                self._mark_excluded(
                    df.index[~positive_mask],
                    f"Non-positive target cannot be transformed with {transform_mode}",
                )
                df = df[positive_mask]
                removed_count = initial_rows - len(df)
                logger.warning(
                    f"Removed {removed_count} rows with non-positive values before {transform_mode} transformation.")

            numeric_vals = pd.to_numeric(df[self.target_col], errors='coerce')
            if transform_mode == "-log10":
                df[self.target_col] = -np.log10(numeric_vals)
            elif transform_mode == "log10":
                df[self.target_col] = np.log10(numeric_vals)
            output_scale = transform_mode

        if self.dropna_conditions and self.condition_cols:
            missing_condition = df[self.condition_cols].isna().any(axis=1)
            self._mark_excluded(
                df.index[missing_condition],
                "Missing required deduplication condition",
            )
            df = df[~missing_condition].copy()

        group_keys = [self.smiles_col] + self.condition_cols
        grouped = df.groupby(group_keys, group_keys=False, sort=False, dropna=False)

        if self.data_type == 'smiles' or not self.target_col:
            result = self._process_without_target(grouped)
        else:
            result = self._process_with_target(grouped, self._p_threshold, progress_callback)
        result.attrs = deepcopy(source_df.attrs)
        if self.target_col and self.data_type == 'continuous':
            # The working frame may be numeric/log-transformed, but source values
            # in the public result must still match the selected original rows.
            result[self.target_col] = source_df[self.target_col].reindex(result.index)
            scales = dict(source_df.attrs.get(TARGET_SCALES_ATTR, {}))
            scales[self.target_col + '_new'] = output_scale
            result.attrs[TARGET_SCALES_ATTR] = scales
        return result

    def _validate_columns(self, df: pd.DataFrame):
        """Validate column existence"""
        required_cols = [self.smiles_col]
        if self.target_col:
            required_cols.append(self.target_col)
        required_cols.extend(self.condition_cols)

        missing = [col for col in required_cols if col not in df.columns]
        if missing:
            raise KeyError(f"Missing required columns: {missing}")

    def _process_without_target(self, grouped) -> pd.DataFrame:
        """Simple deduplication without target value"""
        logger.info("Performing simple deduplication by SMILES")
        processed = []
        for _, group in grouped:
            first_record = group.iloc[[0]].copy()
            processed.append(
                self._mark_record(first_record, method='smiles_only', source_group=group)
            )
        if not processed:
            columns = list(grouped.obj.columns) + [
                'Deduplication Strategy',
                'Deduplication Record Count',
                'Deduplication Source Rows',
                'Deduplication Distinct Value Count',
                'Deduplication Value Range',
            ]
            return pd.DataFrame(columns=list(dict.fromkeys(columns)))
        return pd.concat(processed)

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
                method = f'<=3 -> {self.aggregation}' if self.data_type == 'continuous' else 'no change'
                processed.append(self._mark_record(valid_group, method=method, source_group=valid_group))
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
            new_columns = [
                self.target_col + '_new',
                'Deduplication Strategy',
                'Deduplication Record Count',
                'Deduplication Source Rows',
                'Deduplication Input Values',
                'Deduplication Distinct Value Count',
                'Deduplication Value Range',
            ]
            final_columns = list(dict.fromkeys(original_columns + new_columns))
            return pd.DataFrame(columns=final_columns)

        return pd.concat(processed)

    def _handle_discrete(self, group: pd.DataFrame) -> pd.DataFrame:
        """Handle discrete data"""
        if self.method == 'priority' and self.priority_list:
            for priority_val in self.priority_list:
                matches = group[group[self.target_col] == priority_val]
                if matches.empty:
                    matches = group[group[self.target_col].astype(str) == str(priority_val)]
                if not matches.empty:
                    final_record = matches.head(1).copy()
                    final_record[self.target_col + '_new'] = final_record[self.target_col]
                    return self._mark_record(
                        final_record, method=f'priority_{priority_val}', source_group=group
                    )

        counts = group[self.target_col].value_counts()
        max_count = counts.max()

        if len(counts[counts == max_count]) > 1:
            logger.debug(f"Tie detected in group: {group[self.smiles_col].iloc[0]}")
            self._mark_excluded(group.index, "Unresolved discrete label tie")
            return pd.DataFrame()

        selected = counts.idxmax()
        final_record = group[group[self.target_col] == selected].head(1).copy()
        final_record[self.target_col + '_new'] = selected
        return self._mark_record(final_record, method='vote', source_group=group)

    def _handle_continuous(self, group: pd.DataFrame, p_threshold: float) -> pd.DataFrame:
        """Handle continuous data"""
        n = len(group)
        values = group[self.target_col]

        if self.custom_method:
            clean_values, method = self.custom_method(values)
        elif n <= 3:
            clean_values, method = values, '<=3'
        else:
            clean_values, method = self._remove_outliers(values, p_threshold)

        if clean_values.empty:
            logger.warning(
                f"All values in group for SMILES '{group[self.smiles_col].iloc[0]}' "
                f"were removed as outliers. Skipping this group."
            )
            self._mark_excluded(group.index, "All group values removed as outliers")
            return pd.DataFrame()

        if self.aggregation in {'max', 'min'}:
            final_value = clean_values.max() if self.aggregation == 'max' else clean_values.min()
            selected_index = clean_values[clean_values == final_value].index[0]
            final_record = group.loc[[selected_index]].copy()
            final_record[self.target_col + '_new'] = final_value
            self._mark_excluded(group.index[~group.index.isin(clean_values.index)], f"Outlier removed by {method}")
            return self._mark_record(final_record, method=f'{method} -> {self.aggregation}', source_group=group)

        final_value = self._finite_mean(clean_values)
        self._mark_excluded(group.index[~group.index.isin(clean_values.index)], f"Outlier removed by {method}")
        valid_indices = clean_values.index
        valid_group = group.loc[valid_indices]
        best_idx = self._nearest_index(valid_group[self.target_col], final_value)
        final_record = group.loc[[best_idx]].copy()
        final_record[self.target_col + '_new'] = final_value
        return self._mark_record(final_record, method=f'{method} -> mean', source_group=group)

    @staticmethod
    def _finite_mean(values: pd.Series) -> float:
        """Average finite values without overflowing their intermediate sum."""
        numbers = np.asarray(values, dtype=float)
        if not np.isfinite(numbers).all():
            raise ValueError("Continuous aggregation requires finite numeric values")
        scale = np.max(np.abs(numbers))
        if scale == 0:
            return 0.0
        mean = float(np.mean(numbers / scale) * scale)
        if not np.isfinite(mean):
            raise ValueError("Continuous aggregation produced a non-finite target")
        return mean

    @staticmethod
    def _nearest_index(values: pd.Series, target: float):
        scale = max(float(values.abs().max()), abs(target))
        if scale == 0:
            return values.index[0]
        return (values / scale - target / scale).abs().idxmin()

    @staticmethod
    def _scaled_values(values: pd.Series) -> pd.Series:
        scale = values.abs().max()
        return values / scale if scale else values

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
                stat, p = shapiro(DataDeduplicator._scaled_values(values))
            except ValueError:
                return True
        return p > p_threshold

    @staticmethod
    def _3sigma_filter(values: pd.Series, n_sigma: int = 3) -> pd.Series:
        """3-sigma filtering"""
        scaled = DataDeduplicator._scaled_values(values)
        mean = scaled.mean()
        std = scaled.std()
        lower = mean - n_sigma * std
        upper = mean + n_sigma * std
        return values[(scaled >= lower) & (scaled <= upper)]

    def _mark_excluded(self, indices, reason: str) -> None:
        """Record why source rows cannot contribute a modeling record."""
        for index in indices:
            if index not in self.exclusion_reasons:
                self.exclusion_reasons[index] = reason

    @staticmethod
    def _iqr_filter(values: pd.Series, k: float = 1.5) -> pd.Series:
        """Interquartile range (IQR) filtering"""
        scaled = DataDeduplicator._scaled_values(values)
        q1 = scaled.quantile(0.25)
        q3 = scaled.quantile(0.75)
        iqr = q3 - q1
        lower = q1 - k * iqr
        upper = q3 + k * iqr
        return values[(scaled >= lower) & (scaled <= upper)]

    def _mark_record(self, record: pd.DataFrame,
                     method: Optional[str] = None,
                     source_group: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        """Mark processed record"""
        record = record.copy()
        if method:
            record["Deduplication Strategy"] = method
        if source_group is not None:
            for index in record.index:
                self.source_indices[index] = tuple(source_group.index)
            record["Deduplication Record Count"] = len(source_group)
            record["Deduplication Source Rows"] = ";".join(str(index) for index in source_group.index)
            if self.target_col and self.target_col in source_group.columns:
                values = source_group[self.target_col]
                record["Deduplication Input Values"] = ";".join(str(value) for value in values.tolist())
                record["Deduplication Distinct Value Count"] = values.nunique(dropna=True)
                numeric_values = pd.to_numeric(values, errors='coerce').dropna()
                value_range = (
                    float(numeric_values.max()) - float(numeric_values.min())
                    if not numeric_values.empty else np.nan
                )
                record["Deduplication Value Range"] = (
                    value_range if np.isfinite(value_range) else pd.NA
                )
        return record

    @classmethod
    def create_pipeline(cls, steps: List[Dict[str, Any]]) -> List['DataDeduplicator']:
        """Create a processing pipeline"""
        return [cls(**step_config) for step_config in steps]
