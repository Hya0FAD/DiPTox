"""Shared, typed column-value discovery and row selection."""

import numbers

import pandas as pd


def scalar_value(value):
    """Return a portable scalar; all missing-value representations become None."""
    if not pd.api.types.is_scalar(value):
        raise ValueError("Column filtering requires scalar cell values.")
    if pd.isna(value):
        return None
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, (str, bool, int, float)):
        return value
    if hasattr(value, "isoformat"):
        return value.isoformat()
    return str(value)


def value_key(value):
    value = scalar_value(value)
    # Keep True distinct from 1, and numeric 1 distinct from text "1".
    kind = "number" if isinstance(value, numbers.Number) and not isinstance(value, bool) else type(value).__name__
    return kind, value


def column_values(frame, column):
    if column not in frame.columns:
        raise KeyError(f"Column not found: {column}")
    counts = {}
    for value in frame[column]:
        key = value_key(value)
        if key not in counts:
            counts[key] = {"value": scalar_value(value), "count": 0}
        counts[key]["count"] += 1
    return list(counts.values())


def selection_mask(frame, column, values=None, mode="keep"):
    if column not in frame.columns:
        raise KeyError(f"Column not found: {column}")
    if mode not in {"keep", "remove"}:
        raise ValueError("mode must be 'keep' or 'remove'.")
    if values is not None and not isinstance(values, (list, tuple)):
        raise ValueError("values must be a list or tuple of scalar values, or None.")
    if values is None or len(values) == 0:
        return pd.Series(True, index=frame.index, dtype=bool)
    keys = {value_key(value) for value in values}
    matches = frame[column].map(lambda value: value_key(value) in keys).astype(bool)
    return matches if mode == "keep" else ~matches


def merge_rules(values=None, replacement='other', groups=None):
    """Validate all rules before applying them, matching typed original values."""
    if groups is not None:
        if values is not None and len(values):
            raise ValueError('Use groups or values, not both.')
        if not isinstance(groups, list):
            raise ValueError('groups must be an array of {values, replacement} objects.')
        rules = groups
    else:
        if not isinstance(replacement, str) or not replacement.strip():
            raise ValueError('replacement must be a non-empty text label.')
        if values is not None and not isinstance(values, (list, tuple)):
            raise ValueError('values must be a list or tuple.')
        rules = [{'values': values, 'replacement': replacement}] if values else []
    mapping = {}
    for index, rule in enumerate(rules, 1):
        if not isinstance(rule, dict) or set(rule) != {'values', 'replacement'}:
            raise ValueError(f'Group {index}: expected values and replacement fields.')
        selected, label = rule['values'], rule['replacement']
        if not isinstance(selected, (list, tuple)) or not selected:
            raise ValueError(f'Group {index}: values must be a non-empty array.')
        if not isinstance(label, str) or not label.strip():
            raise ValueError(f'Group {index}: replacement must be non-empty text.')
        for value in selected:
            key = value_key(value)
            if key in mapping and mapping[key] != label:
                raise ValueError(f'Group {index}: value {value!r} belongs to conflicting replacements.')
            mapping[key] = label
    return mapping


def merge_values(frame, column, values=None, replacement='other', output_column=None, groups=None):
    """Create a categorical grouping column without changing source cells or rows."""
    if column not in frame.columns:
        raise KeyError(f'Column not found: {column}')
    mapping = merge_rules(values, replacement, groups)
    if not mapping:
        return frame.copy(), column, 0
    output_column = output_column if output_column is not None else column + ' (Merged)'
    if not isinstance(output_column, str) or not output_column.strip():
        raise ValueError('output_column must be a non-empty column name.')
    if output_column in frame.columns:
        raise ValueError('Output column already exists: ' + output_column)
    result = frame.copy()
    result[output_column] = frame[column].astype(object)
    keys = frame[column].astype(object).map(value_key)
    mask = keys.map(lambda key: key in mapping)
    result.loc[mask, output_column] = keys.loc[mask].map(lambda key: mapping[key])
    generated = list(frame.attrs.get('_diptox_merged_columns', []))
    result.attrs['_diptox_merged_columns'] = generated + [output_column]
    return result, output_column, int(mask.sum())
