"""Preserve imported fields when a processing step produces the same names."""

INITIAL_COLUMNS = {"Canonical SMILES", "Processing Log", "Standardization Status", "Is Valid"}
PREPROCESS_COLUMNS = INITIAL_COLUMNS | {
    "Original Canonical SMILES", "Original Fragment Count", "Final Fragment Count",
    "Original Structure Type", "Final Structure Type", "Is Multi-Component",
    "Original Formal Charge", "Final Formal Charge", "Original Contains Metal",
    "Final Contains Metal", "Original Metal Elements", "Final Metal Elements",
    "Original Molecular Weight", "Standardized Molecular Weight", "Structure Changed",
    "Molecular Weight Changed", "Parent Mapping Count", "Parent Mapping Conflict",
    "Parent Source Structures", "Parent Mapping Status",
}


def source_column_renames(available, source_columns, output_columns):
    """Return deterministic, unused names for source fields about to be replaced."""
    occupied = set(available)
    outputs = set(output_columns)
    renamed = {}
    for column in source_columns:
        if column not in occupied or column not in outputs or column in renamed:
            continue
        candidate = f"{column} (Input)"
        number = 2
        while candidate in occupied or candidate in outputs:
            candidate = f"{column} (Input {number})"
            number += 1
        renamed[column] = candidate
        occupied.add(candidate)
    return renamed
