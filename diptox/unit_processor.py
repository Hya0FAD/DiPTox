# diptox/unit_processor.py
import numpy as np
import pandas as pd
import re
from typing import Optional, Callable, Tuple, Dict
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from .logger import log_manager

logger = log_manager.get_logger(__name__)

TARGET_SCALES_ATTR = '_diptox_target_scales'


def get_target_scale(df: pd.DataFrame, target_col: str) -> str:
    """Return a target's numeric scale; unannotated source columns are linear.

    DataFrame attrs map column names to ``linear``, ``log10`` or ``-log10``.
    Callers creating a new target must copy this map before adding its scale.
    """
    return df.attrs.get(TARGET_SCALES_ATTR, {}).get(target_col, 'linear')


class UnitProcessor:
    """Handles unit conversion using mathematical expressions."""

    def __init__(self, rules: dict = None):
        """
        Initializes the UnitProcessor.
        :param rules: A dictionary of conversion rules. {('mg/L', 'ug/L'): 'x * 1000'}
        """
        self.rules = self._get_default_rules()
        if rules:
            self.rules.update(rules)

    def _get_default_rules(self) -> Dict[Tuple[str, str], str]:
        """
        Defines a set of built-in conversion rules based on common scientific units.
        Derived from pint definitions and common usage.
        """
        rules = {('g/L', 'mg/L'): 'x * 1000', ('mg/L', 'g/L'): 'x / 1000', ('mg/L', 'ug/L'): 'x * 1000', ('mg/L', 'µg/L'): 'x * 1000',
                 ('ug/L', 'mg/L'): 'x / 1000', ('µg/L', 'mg/L'): 'x / 1000', ('mg/L', 'ng/L'): 'x * 1000000',
                 ('ng/L', 'mg/L'): 'x / 1000000', ('ml/L', 'ul/L'): 'x * 1000', ('ul/L', 'ml/L'): 'x / 1000',
                 ('%', 'ppm'): 'x * 10000', ('ppm', '%'): 'x / 10000', ('vol%', 'ppm'): 'x * 10000',
                 ('wt%', 'ppm'): 'x * 10000', ('ppm', 'ppb'): 'x * 1000', ('ppb', 'ppm'): 'x / 1000',
                 ('%', 'g/L'): 'x * 10', ('%', 'mg/L'): 'x * 10000', ('d', 'h'): 'x * 24', ('day', 'h'): 'x * 24',
                 ('h', 'd'): 'x / 24', ('h', 'min'): 'x * 60', ('min', 'h'): 'x / 60', ('d', 'min'): 'x * 1440',
                 ('min', 's'): 'x * 60', ('s', 'min'): 'x / 60', ('mmHg', 'Pa'): 'x * 133.322',
                 ('mm Hg', 'Pa'): 'x * 133.322', ('Torr', 'Pa'): 'x * 133.322', ('psi', 'Pa'): 'x * 6894.76',
                 ('atm', 'Pa'): 'x * 101325', ('kPa', 'Pa'): 'x * 1000', ('°C', 'K'): 'x + 273.15',
                 ('degC', 'K'): 'x + 273.15', ('K', '°C'): 'x - 273.15', ('°F', '°C'): '(x - 32) * 5/9',
                 ('degF', '°C'): '(x - 32) * 5/9', ('°C', '°F'): '(x * 9/5) + 32', ('mM', 'M'): 'x / 1000',
                 ('uM', 'M'): 'x / 1000000', ('µM', 'M'): 'x / 1000000', ('nM', 'M'): 'x / 1000000000',
                 ('M', 'mM'): 'x * 1000', ('M', 'uM'): 'x * 1000000', ('M', 'g/L'): 'x * mw',
                 ('M', 'mg/L'): 'x * mw * 1000', ('mM', 'g/L'): '(x * mw) / 1000', ('mM', 'mg/L'): 'x * mw',
                 ('uM', 'mg/L'): '(x * mw) / 1000', ('uM', 'ug/L'): 'x * mw', ('µM', 'ug/L'): 'x * mw',
                 ('nM', 'ng/L'): 'x * mw', ('g/L', 'M'): 'x / mw', ('mg/L', 'M'): 'x / (mw * 1000)',
                 ('g/L', 'mM'): '(x * 1000) / mw', ('mg/L', 'mM'): 'x / mw', ('mg/L', 'uM'): '(x * 1000) / mw',
                 ('ug/L', 'uM'): 'x / mw', ('µg/L', 'µM'): 'x / mw', ('ng/L', 'nM'): 'x / mw'}
        return rules

    def add_rule(self, from_unit: str, to_unit: str, formula: str):
        """Adds or updates a conversion rule."""
        if not self._is_valid_formula(formula):
            raise ValueError(f"Invalid or unsafe formula provided: {formula}")
        self.rules[(from_unit, to_unit)] = formula
        logger.info(f"Added conversion rule: {from_unit} -> {to_unit} | {formula}")

    def get_rule(self, from_unit: str, to_unit: str) -> str or None:
        """Retrieves a conversion rule."""
        return self.rules.get((from_unit, to_unit))

    @staticmethod
    def _is_valid_formula(formula: str) -> bool:
        """
        Validates the formula to ensure it only contains allowed elements.
        - 'x' as the variable
        - Numbers (integers and floats, including scientific notation)
        - Basic operators: +, -, *, /, **
        - Parentheses: ()
        - Allowed functions: log, log10, exp
        """
        if 'x' not in formula:
            logger.error(f"Formula must contain 'x' as the variable: {formula}")
            return False

        # Allow only specific characters and patterns
        allowed_chars_pattern = r"^[x\d\s\.\+\-\*\/\(\)e\^]+$"
        clean_formula = formula.replace("log10", "").replace("log", "").replace("exp", "").replace("mw", "")
        if not re.match(allowed_chars_pattern, clean_formula):
            logger.error(f"Formula contains disallowed characters: {formula}")
            return False

        # Prevent calling other functions or using other variables
        # Finds any word that is not 'x', 'mw', 'log', 'log10', 'exp', or 'e'
        disallowed_names = re.findall(r"\b(?!x|mw|log10|log|exp|e\b)[a-df-zA-Z_]\w*\b", formula)
        if disallowed_names:
            logger.error(f"Formula contains disallowed names: {disallowed_names}")
            return False

        return True

    @staticmethod
    def _apply_precision(original_series: pd.Series, converted_series: pd.Series) -> pd.Series:
        def get_sig_figs(s):
            if pd.isna(s):
                return np.nan
            s = str(s).strip().lower().lstrip('-')
            if 'e' in s:
                s = s.split('e')[0]
            s_clean = s.replace('.', '').lstrip('0')
            if not s_clean:
                return 1
            return len(s_clean)

        def round_to_sig_figs(val, sf):
            if pd.isna(val) or pd.isna(sf):
                return val
            if val == 0:
                return 0.0
            sf = int(sf)
            try:
                order = int(np.floor(np.log10(abs(val))))
                decimals = sf - 1 - order
                if decimals <= 0:
                    return np.round(val, 0)
                return np.round(val, decimals)
            except Exception:
                return val
        sig_figs = original_series.apply(get_sig_figs)
        rounded_values = [round_to_sig_figs(v, sf) for v, sf in zip(converted_series, sig_figs)]
        return pd.Series(rounded_values, index=converted_series.index)

    def convert(self, values: pd.Series, formula: str, mw: Optional[pd.Series] = None) -> pd.Series:
        """
        Applies the conversion formula to a pandas Series.
        :param values: The series of numerical data to convert.
        :param formula: The mathematical expression for conversion, using 'x' as the variable.
        :param mw: An optional series of molecular weights, aligned with 'values'.
        :return: A new Series with the converted values.
        """
        if not self._is_valid_formula(formula):
            raise ValueError(f"Invalid or unsafe formula provided for conversion: {formula}")

        python_formula = formula.replace('^', '**')

        # pandas 2.0 (Python 3.8) cannot evaluate arithmetic with nullable Float64.
        # Normalize calculation inputs to NumPy floats while preserving missing values.
        numeric_values = pd.to_numeric(values, errors='coerce').astype(float)
        numeric_values = numeric_values.where(np.isfinite(numeric_values))

        safe_dict = {
            'x': numeric_values,
            'log': np.log,
            'log10': np.log10,
            'exp': np.exp,
            'e': np.e,
        }

        if 'mw' in formula:
            if mw is None:
                raise ValueError("Formula requires 'mw' but no molecular weight series was provided.")
            numeric_mw = pd.to_numeric(mw, errors='coerce').astype(float)
            safe_dict['mw'] = numeric_mw.where(np.isfinite(numeric_mw) & numeric_mw.gt(0))

        try:
            result = pd.eval(python_formula, local_dict=safe_dict, global_dict={})
            result = pd.to_numeric(pd.Series(result, index=values.index), errors='coerce')
            usable = np.isfinite(result) & numeric_values.notna()
            if 'mw' in formula:
                usable &= safe_dict['mw'].notna()
            return result.where(usable).astype(float)
        except Exception as e:
            logger.error(f"Failed to evaluate formula '{formula}': {e}")
            return pd.Series(np.nan, index=values.index)

    def standardize(self, df: pd.DataFrame, target_col: str, unit_col: str, standard_unit: str,
                    smiles_col: Optional[str] = None,
                    molecular_weight_col: Optional[str] = None,
                    rule_provider_callback: Optional[Callable[[str, str], Optional[str]]] = None) -> Tuple[pd.DataFrame, str, str]:
        """
        Performs the full unit standardization process on a DataFrame.

        :param df: The DataFrame to process.
        :param target_col: The name of the column with values to convert.
        :param unit_col: The name of the column specifying the units.
        :param standard_unit: The target unit to convert all values to.
        :param smiles_col: The name of the column containing SMILES strings, for MW-based conversions.
        :param molecular_weight_col: Optional explicit molecular-weight column. This takes precedence over SMILES.
        :param rule_provider_callback: An optional function that takes (from_unit, to_unit)
                                       and returns a formula string if a rule is missing.
        :return: A tuple containing the processed DataFrame and the new target column name.
        """
        scale = get_target_scale(df, target_col)
        source_units = df[unit_col].dropna().astype(str)
        has_log_units = source_units.str.match(r'^-?log10\(.*\)$').any()
        if scale != 'linear' or has_log_units:
            raise ValueError(
                f"Target '{target_col}' is on a logarithmic scale ({scale}); "
                "linear unit conversion is not defined. Use the original linear "
                "target or explicitly inverse-transform it first."
            )
        df = df.copy()
        numeric_values = pd.to_numeric(df[target_col], errors='coerce')
        invalid_mask = numeric_values.isna() | ~np.isfinite(numeric_values)

        unique_units = [u for u in df[unit_col].dropna().unique() if u]
        new_target_col = f"{target_col} (Standardized)"
        new_unit_col = f"{unit_col} (Standardized)"
        status_col = "Unit Conversion Status"
        mw_used_col = "Molecular Weight Used"
        mw_source_col = "Molecular Weight Source"
        df[new_target_col] = pd.NA
        df[new_unit_col] = pd.NA
        df[status_col] = pd.Series(pd.NA, index=df.index, dtype="string")
        df[mw_used_col] = pd.Series(pd.NA, index=df.index, dtype="Float64")
        df[mw_source_col] = pd.Series(pd.NA, index=df.index, dtype="string")
        df.loc[invalid_mask, status_col] = "Invalid numeric value"
        df.loc[df[unit_col].isna() | df[unit_col].eq(""), status_col] = "Missing source unit"

        if not standard_unit:
            raise ValueError("A standard unit must be provided.")
        if standard_unit not in unique_units:
            logger.info(f"Converting all source units to the new standard unit '{standard_unit}'.")

        for unit in unique_units:
            mask = (df[unit_col] == unit) & ~invalid_mask
            if not mask.any():
                continue
            if unit == standard_unit:
                df.loc[mask, new_target_col] = numeric_values.loc[mask]
                df.loc[mask, status_col] = "No conversion needed"
                df.loc[mask, mw_source_col] = "not required"
                continue

            formula = self.get_rule(unit, standard_unit)
            if not formula and rule_provider_callback:
                formula_input = rule_provider_callback(unit, standard_unit)
                if formula_input:
                    try:
                        self.add_rule(unit, standard_unit, formula_input)
                        formula = formula_input
                    except ValueError as e:
                        logger.error(f"Invalid formula provided: {e}. Skipping unit '{unit}'.")
                        formula = None
                else:
                    logger.warning(f"Skipping conversion for unit '{unit}' as no rule was provided.")

            if formula:
                values_to_convert = numeric_values.loc[mask]
                mw_values = None
                if 'mw' in formula:
                    if molecular_weight_col:
                        if molecular_weight_col not in df.columns:
                            raise ValueError(
                                f"Molecular-weight column '{molecular_weight_col}' was not found."
                            )
                        mw_values = pd.to_numeric(
                            df.loc[mask, molecular_weight_col], errors='coerce'
                        )
                        mw_source = f"column:{molecular_weight_col}"
                    else:
                        if not smiles_col or smiles_col not in df.columns:
                            raise ValueError(
                                "Formula requires an explicit molecular-weight basis. Choose "
                                "molecular_weight_source='original' or 'standardized' in the pipeline, "
                                "or provide molecular_weight_col (or an explicit smiles_col directly)."
                            )

                        def calculate_mw(smiles):
                            if not isinstance(smiles, str) or not smiles.strip():
                                return None
                            try:
                                mol = Chem.MolFromSmiles(smiles)
                                return Descriptors.MolWt(mol) if mol else None
                            except Exception:
                                return None

                        mw_values = df.loc[mask, smiles_col].apply(calculate_mw)
                        mw_source = f"structure:{smiles_col}"

                    mw_values = pd.to_numeric(mw_values, errors='coerce')
                    available_mw = (mw_values.notna() & np.isfinite(mw_values) & mw_values.gt(0)).fillna(False)
                    df.loc[mw_values.index[~available_mw], status_col] = "Invalid molecular weight"
                    mw_values = mw_values.where(available_mw)
                    df.loc[mw_values.index[available_mw], mw_used_col] = mw_values[available_mw]
                    df.loc[mw_values.index[available_mw], mw_source_col] = mw_source
                else:
                    df.loc[mask, mw_source_col] = "not required"

                converted_values = self.convert(values_to_convert, formula, mw=mw_values)
                converted_values = self._apply_precision(df.loc[mask, target_col], converted_values)
                successful = converted_values.notna() & np.isfinite(converted_values)
                successful_indices = converted_values.index[successful]
                df.loc[successful_indices, new_target_col] = converted_values[successful]
                df.loc[successful_indices, status_col] = "Converted"
                failed_indices = converted_values.index[~successful]
                unset_status = df.loc[failed_indices, status_col].isna()
                df.loc[unset_status[unset_status].index, status_col] = "Conversion failed"
            else:
                df.loc[mask, status_col] = "Missing conversion rule"
                if not rule_provider_callback:
                    logger.warning(f"Missing conversion rule for '{unit}' -> '{standard_unit}'.")

        # Finalize
        df.loc[df[new_target_col].notna(), new_unit_col] = standard_unit
        scales = dict(df.attrs.get(TARGET_SCALES_ATTR, {}))
        scales[new_target_col] = 'linear'
        df.attrs[TARGET_SCALES_ATTR] = scales
        return df, new_target_col, new_unit_col
