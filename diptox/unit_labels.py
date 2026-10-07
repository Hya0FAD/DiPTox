"""Display logarithmic units as text when opening CSV/TXT files in Excel."""

import re


def normalize_unit_label(value):
    if not isinstance(value, str):
        return value
    match = re.fullmatch(r"=?([-\N{MINUS SIGN}]?log10)\((.*)\)", value, re.IGNORECASE)
    if match:
        return match.group(1).lower().replace('\N{MINUS SIGN}', '-') + '(' + match.group(2) + ')'
    return value


def excel_unit_label(value):
    normalized = normalize_unit_label(value)
    if isinstance(normalized, str) and normalized.startswith('-log10('):
        return '\N{MINUS SIGN}' + normalized[1:]
    return normalized
