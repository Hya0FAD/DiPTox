"""Unit and base-10 transformations for independent measurement columns."""
import numpy as np
import pandas as pd
from copy import deepcopy

from .unit_labels import normalize_unit_label
from .unit_processor import UnitProcessor, TARGET_SCALES_ATTR, get_target_scale


def output_columns(value_col, unit_col, log_transform):
    suffix = ' (Transformed)' if log_transform != 'None' else ' (Standardized)'
    return value_col + suffix, unit_col + suffix


def transform_frame(frame, value_col, unit_col, standard_unit=None,
                    log_transform='None', conversion_rules=None,
                    smiles_col=None, molecular_weight_col=None):
    if log_transform not in ('None', 'log10', '-log10'):
        raise ValueError('log_transform must be None (string), log10 or -log10.')
    if value_col == unit_col:
        raise ValueError('Select distinct value and unit columns.')
    for column in (value_col, unit_col, molecular_weight_col):
        if column is not None and column not in frame:
            raise ValueError('Column not found: ' + column)
    if not standard_unit and log_transform == 'None':
        raise ValueError('Specify a standard_unit or a log_transform.')
    new_value, new_unit = output_columns(value_col, unit_col, log_transform)
    status = new_value + ' Status'
    for column in (new_value, new_unit, status):
        if column in frame:
            raise ValueError('Output column already exists; undo or select another source: ' + column)
    units = frame[unit_col].map(normalize_unit_label)
    if (get_target_scale(frame, value_col) != 'linear'
            or units.dropna().astype(str).str.match(r'^-?log10\(').any()):
        raise ValueError('Select original linear values before unit/log transformation.')
    if standard_unit and str(standard_unit).startswith(('log10(', '-log10(', '−log10(')):
        raise ValueError('Use log_transform to select a logarithmic scale.')
    processor = UnitProcessor(conversion_rules)
    work = frame[[value_col, unit_col]].copy()
    work[unit_col] = units
    for column in (smiles_col, molecular_weight_col):
        if column and column in frame:
            work[column] = frame[column]
    if standard_unit:
        for unit in units.dropna().unique():
            if not unit or unit == standard_unit:
                continue
            rule = processor.get_rule(unit, standard_unit)
            if not rule:
                raise ValueError('Missing conversion rule: %s -> %s' % (unit, standard_unit))
            if not processor._is_valid_formula(rule):
                raise ValueError('Invalid conversion formula: ' + rule)
            if 'mw' in rule and not (smiles_col or molecular_weight_col):
                raise ValueError('This conversion requires an explicit molecular-weight basis or column.')
        work, converted_value, converted_unit = processor.standardize(
            work, value_col, unit_col, standard_unit,
            smiles_col=smiles_col, molecular_weight_col=molecular_weight_col)
        values = pd.to_numeric(work[converted_value], errors='coerce').astype(float)
        result_units = work[converted_unit]
        reasons = work['Unit Conversion Status'].copy()
    else:
        values = pd.to_numeric(frame[value_col], errors='coerce').astype(float)
        result_units = units.copy()
        reasons = pd.Series('No conversion needed', index=frame.index)
    valid = np.isfinite(values) & result_units.notna() & result_units.ne('')
    reasons.loc[~valid & reasons.isin(['Converted', 'No conversion needed'])] = 'Invalid value or missing unit'
    if log_transform != 'None':
        nonpositive = valid & values.le(0)
        reasons.loc[nonpositive] = 'Log transformation requires a positive value'
        valid &= values.gt(0)
        values = np.log10(values.where(valid))
        if log_transform == '-log10':
            values = -values
        result_units = result_units.map(lambda unit: '%s(%s)' % (log_transform, unit))
    result = frame.copy()
    result[new_value] = values.where(valid)
    result[new_unit] = result_units.where(valid)
    result[status] = reasons
    scales = dict(frame.attrs.get(TARGET_SCALES_ATTR, {}))
    scales[new_value] = 'linear' if log_transform == 'None' else log_transform
    result.attrs[TARGET_SCALES_ATTR] = scales
    dependencies = deepcopy(frame.attrs.get('_diptox_standardized_mw_dependencies', {}))
    inherited = dict(dependencies.get(value_col, {}))
    if standard_unit and molecular_weight_col == 'Standardized Molecular Weight':
        for index, weight in work['Molecular Weight Used'].items():
            if pd.notna(weight):
                inherited[index] = float(weight)
    if inherited:
        dependencies[new_value] = inherited
    result.attrs['_diptox_standardized_mw_dependencies'] = dependencies
    pairs = dict(frame.attrs.get('_diptox_condition_columns', {}))
    pairs[new_value] = new_unit
    result.attrs['_diptox_condition_columns'] = pairs
    return result, new_value, new_unit, valid, reasons
