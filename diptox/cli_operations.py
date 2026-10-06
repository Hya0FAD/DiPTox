"""Offline CLI operations and per-run chemistry rules."""

from copy import deepcopy

import pandas as pd
from rdkit import Chem

from .cli_config import CliError


def _rule_error(field, message, value=None):
    raise CliError("INVALID_RULE", message, {"field": field, "value": value})


def _rule_molecule(value, field, smarts=False, sanitize=True):
    if not isinstance(value, str) or not value.strip():
        _rule_error(field, "A nonempty molecular pattern is required.", value)
    try:
        mol = Chem.MolFromSmarts(value) if smarts else Chem.MolFromSmiles(value, sanitize=sanitize)
        if mol is not None and not smarts and not sanitize:
            # Replacement fragments may contain non-ring aromatic atoms, but
            # impossible valences should still fail during configuration.
            mol.UpdatePropertyCache(strict=True)
            Chem.FastFindRings(mol)
    except Exception:
        mol = None
    if mol is None or not mol.GetNumAtoms():
        _rule_error(field, "Invalid SMARTS pattern." if smarts else "Invalid SMILES pattern.", value)
    return mol


def apply_rules(processor, rules):
    """Apply validated rule changes to one processor, without persistent state.

    Salts are SMARTS, solvents are SMILES. Neutralization replacement fragments
    are parsed without sanitization, as in ChemistryProcessor (e.g. aromatic n).
    A failed change leaves the processor untouched. Removals precede additions.
    """
    if not isinstance(rules, dict) or set(rules) - {"atoms", "salts", "solvents", "neutralization"}:
        _rule_error("rules", "Unknown or invalid rule sections.")
    candidate = deepcopy(processor)
    symbols = {Chem.GetPeriodicTable().GetElementSymbol(number) for number in range(1, 119)}
    for section, changes in rules.items():
        if not isinstance(changes, dict) or set(changes) - {"add", "remove"}:
            _rule_error(f"rules.{section}", "Rule changes must contain only add and remove lists.")
        for action in ("remove", "add"):
            values = changes.get(action, [])
            if not isinstance(values, list):
                _rule_error(f"rules.{section}.{action}", "Rule changes must be lists.")
            for number, value in enumerate(values):
                field = f"rules.{section}.{action}[{number}]"
                if section == "atoms":
                    if not isinstance(value, str) or value not in symbols:
                        _rule_error(field, "Invalid element symbol.", value)
                    if action == "remove":
                        if value not in candidate._valid:
                            _rule_error(field, "Cannot remove an atom absent from the effective rules.", value)
                        candidate._valid.remove(value)
                    else:
                        candidate._valid.add(value)
                elif section == "neutralization":
                    if action == "remove":
                        _rule_molecule(value, field, smarts=True)
                        if not candidate.remove_neutralization_rule(value):
                            _rule_error(field, "Cannot remove an absent neutralization rule.", value)
                    else:
                        if not isinstance(value, dict) or set(value) != {"reactant", "product"}:
                            _rule_error(field, "A neutralization rule requires reactant and product.", value)
                        reactant = _rule_molecule(value["reactant"], f"{field}.reactant", smarts=True)
                        product = _rule_molecule(value["product"], f"{field}.product", sanitize=False)
                        if product.HasSubstructMatch(reactant):
                            _rule_error(field, "Neutralization replacement still matches its reactant and may never terminate.", value)
                        if not candidate.add_neutralization_rule(value["reactant"], value["product"]):
                            _rule_error(field, "Could not apply the neutralization rule.", value)
                else:
                    mol = _rule_molecule(value, field, smarts=section == "salts")
                    canonical = Chem.MolToSmarts if section == "salts" else Chem.MolToSmiles
                    key = canonical(mol)
                    effective = getattr(candidate, f"_get_effective_{section}")()
                    removed = getattr(candidate, f"_removed_{section}")
                    custom = getattr(candidate, f"_custom_{section}")
                    if action == "remove":
                        if key not in {canonical(item) for item in effective}:
                            _rule_error(field, "Cannot remove a pattern absent from the effective rules.", value)
                        removed.append(mol)
                    else:
                        removed[:] = [item for item in removed if canonical(item) != key]
                        if key not in {canonical(item) for item in effective}:
                            custom.append(mol)
    for attribute in ("_valid", "_neutralization_rules", "_custom_salts", "_removed_salts",
                      "_custom_solvents", "_removed_solvents"):
        setattr(processor, attribute, getattr(candidate, attribute))
    current = processor.get_current_rules_dict()
    current["neutralization"] = [
        {"reactant": reactant, "product": product}
        for reactant, product in current["neutralization"]
    ]
    return current


def preflight_operation(op, params):
    """Validate chemical query syntax without processing input molecules."""
    if op != "search":
        return
    pattern = params.get("query_pattern")
    query_type = params.get("query_type", "smarts")
    try:
        mol = None
        if isinstance(pattern, str) and pattern.strip():
            mol = Chem.MolFromSmarts(pattern) if query_type == "smarts" else Chem.MolFromSmiles(pattern)
    except Exception:
        mol = None
    if mol is None or not mol.GetNumAtoms():
        raise CliError("INVALID_QUERY", f"Invalid {query_type.upper()} query pattern.",
                       {"query_pattern": pattern, "query_type": query_type})


def _molecules(pipeline, before):
    column = "Canonical SMILES" if pipeline._preprocess_key else pipeline.smiles_col
    molecules, errors = [], []
    for _, row in before.iterrows():
        value = row[column]
        reason = None
        if not isinstance(value, str):
            missing = pd.api.types.is_scalar(value) and pd.isna(value)
            reason = "Missing or empty SMILES" if missing else "Invalid SMILES"
        elif not value.strip():
            reason = "Missing or empty SMILES"
        elif pipeline._preprocess_key and pd.notna(row.get("Is Valid")) and not row["Is Valid"]:
            reason = "Structure rejected by preprocessing"
        mol = None if reason else pipeline.chem_processor.smiles_to_mol(value, sanitize=True)
        if reason is None and (mol is None or not mol.GetNumAtoms()):
            reason = "Invalid SMILES"
        molecules.append(mol if reason is None else None)
        errors.append(reason)
    return column, molecules, errors


def _audit(before, reasons, severities):
    mask = [reason is not None for reason in reasons]
    audit = before.loc[mask].copy()
    audit["CLI Reason"] = [reason for reason in reasons if reason is not None]
    audit["CLI Severity"] = [severity for reason, severity in zip(reasons, severities) if reason is not None]
    return audit


def _search(pipeline, params, before, column, errors):
    pattern = params["query_pattern"]
    result_column = f"Substructure_{pattern}"
    pipeline.df[result_column] = pd.Series(False, index=pipeline.df.index, dtype="boolean")
    # The historical searcher directly feeds raw values to RDKit. Skip invalid
    # inputs there and report them explicitly as unknown, never as nonmatches.
    source = pipeline.df[column].copy()
    pipeline.df[column] = pd.Series(
        [str(value) if error is None else pd.NA for value, error in zip(source, errors)],
        index=pipeline.df.index, dtype="string")
    try:
        pipeline.substructure_search(pattern, is_smarts=params.get("query_type", "smarts") == "smarts")
    finally:
        pipeline.df[column] = source
    invalid = pd.Series([reason is not None for reason in errors], index=before.index, dtype=bool)
    pipeline.df.loc[invalid, result_column] = pd.NA
    matches = pipeline.df[result_column].fillna(False)
    mode = params.get("mode", "annotate")
    reasons = list(errors)
    severities = ["error" if reason else None for reason in errors]
    filtered = 0
    if mode != "annotate":
        selected = (matches if mode == "matches" else ~matches) & ~invalid
        for position, keep in enumerate(selected):
            if not keep and reasons[position] is None:
                reasons[position] = "Substructure does not match query" if mode == "matches" else "Substructure matches query"
                severities[position] = "filtered"
                filtered += 1
        pipeline.df = pipeline.df.loc[selected].copy()
        pipeline._record_step("Search Selection", before, pipeline.df, f"Mode: {mode} | Pattern: {pattern}")
    return _audit(before, reasons, severities), {
        "matched_rows": int(matches.sum()), "nonmatched_rows": int((~matches & ~invalid).sum()),
        "invalid_rows": int(invalid.sum()), "rejected_rows": int(invalid.sum()), "filtered_rows": filtered,
    }


def _filter_atoms(pipeline, params, before, molecules, errors):
    reasons = list(errors)
    severities = ["error" if reason else None for reason in errors]
    for position, mol in enumerate(molecules):
        if mol is None:
            continue
        counts = {"heavy": mol.GetNumHeavyAtoms()}
        if params.get("min_total_atoms") is not None or params.get("max_total_atoms") is not None:
            try:
                counts["total"] = Chem.AddHs(mol).GetNumAtoms()
            except Exception:
                reasons[position] = "Unable to calculate total atom count"
                severities[position] = "error"
                continue
        failures = []
        for kind, count in counts.items():
            minimum, maximum = params.get(f"min_{kind}_atoms"), params.get(f"max_{kind}_atoms")
            if minimum is not None and count < minimum:
                failures.append(f"{kind.capitalize()} atom count {count} is below minimum {minimum}")
            if maximum is not None and count > maximum:
                failures.append(f"{kind.capitalize()} atom count {count} exceeds maximum {maximum}")
        if failures:
            reasons[position] = "; ".join(failures)
            severities[position] = "filtered"
    if before.empty:
        pipeline._record_step("Filter Atom Count", before, pipeline.df, "No input rows")
    else:
        pipeline.filter_by_atom_count(**params)
    failures = severities.count("error")
    return _audit(before, reasons, severities), {
        "invalid_rows": failures, "rejected_rows": failures, "filtered_rows": severities.count("filtered"),
    }


def _inchi(pipeline, before, column, errors):
    pipeline.df["InChI"] = pd.Series(pd.NA, index=pipeline.df.index, dtype="string")
    source = pipeline.df[column].copy()
    pipeline.df[column] = pd.Series(
        [value if error is None else pd.NA for value, error in zip(source, errors)],
        index=pipeline.df.index, dtype="string")
    try:
        pipeline.calculate_inchi()
    finally:
        pipeline.df[column] = source
    generated = pipeline.df["InChI"].fillna("").str.startswith("InChI=")
    reasons = [reason or (None if success else "InChI generation failed")
               for reason, success in zip(errors, generated)]
    failures = sum(reason is not None for reason in reasons)
    return _audit(before, reasons, ["error"] * len(reasons)), {
        "generated_rows": int(generated.sum()), "invalid_rows": sum(reason is not None for reason in errors),
        "rejected_rows": failures, "filtered_rows": 0,
    }


def execute_operation(pipeline, op, params):
    """Run an offline extension, returning source-row audit records and counts."""
    preflight_operation(op, params)
    source_column = "Canonical SMILES" if pipeline._preprocess_key else pipeline.smiles_col
    result_column = "InChI" if op == "inchi" else f"Substructure_{params['query_pattern']}" if op == "search" else None
    if result_column is not None and result_column == source_column:
        raise CliError("RESULT_COLUMN_CONFLICT", "The generated result column would overwrite the active SMILES column.",
                       {"column": source_column, "operation": op})
    if result_column is not None:
        pipeline._protect_source_columns({result_column})
    before = pipeline.df.copy()
    column, molecules, errors = _molecules(pipeline, before)
    if op == "search":
        return _search(pipeline, params, before, column, errors)
    if op == "filter-atoms":
        return _filter_atoms(pipeline, params, before, molecules, errors)
    if op == "inchi":
        return _inchi(pipeline, before, column, errors)
    raise CliError("INVALID_CONFIG", f"Unsupported offline operation: {op}")
