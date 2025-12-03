# diptox/unit_processor.py
import numpy as np
import pandas as pd
import re
from typing import Optional, Callable, Tuple, Dict
from .logger import log_manager

logger = log_manager.get_logger(__name__)


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
        rules = {('g/L', 'mg/L'): 'x * 1000', ('mg/L', 'g/L'): 'x / 1000', ('mg/L', 'ug/L'): 'x * 1000',
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
                 ('M', 'mM'): 'x * 1000', ('M', 'uM'): 'x * 1000000'}
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
        allowed_chars_pattern = r"^[x\d\s\.\+\-\*\/\(\)e]+$"
        if not re.match(allowed_chars_pattern, formula.replace("log10", "").replace("log", "").replace("exp", "")):
            logger.error(f"Formula contains disallowed characters: {formula}")
            return False

        # Prevent calling other functions or using other variables
        # Finds any word that is not 'x', 'log', 'log10', 'exp', or 'e'
        disallowed_names = re.findall(r"\b(?!x|log10|log|exp|e\b)[a-df-zA-Z_]\w*\b", formula)
        if disallowed_names:
            logger.error(f"Formula contains disallowed names: {disallowed_names}")
            return False

        return True

    def convert(self, values: pd.Series, formula: str) -> pd.Series:
        """
        Applies the conversion formula to a pandas Series.
        :param values: The series of numerical data to convert.
        :param formula: The mathematical expression for conversion, using 'x' as the variable.
        :return: A new Series with the converted values.
        """
        if not self._is_valid_formula(formula):
            raise ValueError(f"Invalid or unsafe formula provided for conversion: {formula}")

        numeric_values = pd.to_numeric(values, errors='coerce')

        safe_dict = {
            'x': numeric_values,
            'log': np.log,
            'log10': np.log10,
            'exp': np.exp,
            'e': np.e,
        }

        try:
            result = pd.eval(formula, local_dict=safe_dict, global_dict={})
            return result
        except Exception as e:
            logger.error(f"Failed to evaluate formula '{formula}': {e}")
            return pd.Series(np.nan, index=values.index)

    def standardize(self, df: pd.DataFrame, target_col: str, unit_col: str, standard_unit: str,
                    rule_provider_callback: Optional[Callable[[str, str], Optional[str]]] = None) -> Tuple[pd.DataFrame, str, str]:
        """
        Performs the full unit standardization process on a DataFrame.

        :param df: The DataFrame to process.
        :param target_col: The name of the column with values to convert.
        :param unit_col: The name of the column specifying the units.
        :param standard_unit: The target unit to convert all values to.
        :param rule_provider_callback: An optional function that takes (from_unit, to_unit)
                                       and returns a formula string if a rule is missing.
        :return: A tuple containing the processed DataFrame and the new target column name.
        """
        unique_units = [u for u in df[unit_col].dropna().unique() if u]
        new_target_col = f"{target_col} (Standardized)"
        new_unit_col = f"{unit_col} (Standardized)"
        df[new_target_col] = pd.NA
        df[new_unit_col] = pd.NA

        if not standard_unit:
            raise ValueError("A standard unit must be provided.")
        if standard_unit not in unique_units:
            logger.warning(f"The specified standard unit '{standard_unit}' is not present in the data.")

        for unit in unique_units:
            mask = df[unit_col] == unit
            if unit == standard_unit:
                df.loc[mask, new_target_col] = df.loc[mask, target_col]
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
                values_to_convert = df.loc[mask, target_col]
                converted_values = self.convert(values_to_convert, formula)
                df.loc[mask, new_target_col] = converted_values
            else:
                if not rule_provider_callback:
                    logger.warning(f"Missing conversion rule for '{unit}' -> '{standard_unit}'. Data for this unit will be removed.")

        # Finalize
        initial_rows = len(df)
        df.loc[df[new_target_col].notna(), new_unit_col] = standard_unit
        df.dropna(subset=[new_target_col], inplace=True)
        final_rows = len(df)
        if initial_rows > final_rows:
            logger.warning(f"{initial_rows - final_rows} rows removed due to failed or skipped unit conversions.")

        return df, new_target_col, new_unit_col
