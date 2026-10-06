"""File-based, non-interactive execution shared by all processing commands."""

from __future__ import annotations

import ast
import hashlib
import json
import math
import os
import platform
import re
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

from .cli_config import CliError
from .column_policy import INITIAL_COLUMNS, PREPROCESS_COLUMNS, source_column_renames


SOURCE_ROWS = "Diptox Source Rows"
AUDIT_COLUMNS = [SOURCE_ROWS, "CLI Step", "CLI Operation", "CLI Reason", "CLI Severity"]


def envelope(command, **values):
    from . import __version__
    return {
        "schema_version": "1", "diptox_version": __version__,
        "status": "success", "command": command, "summary": {},
        "artifacts": {}, "warnings": [], "error": None, **values,
    }


def json_safe(value):
    """Normalize pandas/NumPy scalars without emitting non-standard NaN JSON."""
    import pandas as pd
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if value is None or isinstance(value, (str, bool)):
        return value
    if pd.isna(value):
        return None
    if hasattr(value, "item"):
        return json_safe(value.item())
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, int):
        return value
    if hasattr(value, "isoformat"):
        return value.isoformat()
    return str(value)


def _sheets(path):
    import pandas as pd
    with pd.ExcelFile(path) as book:
        return book.sheet_names


def read_frame(spec):
    """Read the explicitly selected table, never prompting for a sheet."""
    import pandas as pd
    path = Path(spec["path"])
    if not path.is_file():
        raise CliError("INPUT_NOT_FOUND", "Input file does not exist.", {"path": str(path)}, exit_code=3)
    suffix = path.suffix.lower()
    header = 0 if spec.get("header", True) else None
    if "sheet" in spec and suffix not in {".xlsx", ".xls"}:
        raise CliError("INVALID_ARGUMENTS", "A worksheet can only be selected for Excel input.")
    if "delimiter" in spec and suffix not in {".csv", ".txt", ".smi"}:
        raise CliError("INVALID_ARGUMENTS", "A delimiter can only be selected for text input.")
    if suffix in {".xlsx", ".xls"}:
        sheets = _sheets(path)
        sheet = spec.get("sheet")
        if sheet is None:
            if len(sheets) != 1:
                raise CliError("SHEET_REQUIRED", "Select a worksheet using --sheet or input.sheet.", {"sheets": sheets})
            sheet = sheets[0]
        if (isinstance(sheet, int) and not 0 <= sheet < len(sheets)) or (isinstance(sheet, str) and sheet not in sheets):
            raise CliError("INVALID_SHEET", "Worksheet was not found.", {"sheet": sheet, "sheets": sheets})
        frame = pd.read_excel(path, sheet_name=sheet, header=header)
    elif suffix in {".csv", ".txt"}:
        frame = pd.read_csv(path, header=header, encoding=spec.get("encoding", "utf-8"),
                            sep=spec.get("delimiter", "\t" if suffix == ".txt" else ","))
    elif suffix == ".smi":
        from .data_io import DataHandler
        columns = spec.get("columns", {})
        options = {"sep": spec["delimiter"]} if "delimiter" in spec else {}
        frame = DataHandler._load_smi(str(path), smiles_col=columns.get("smiles", "smiles"),
                                     id_col=columns.get("id"), header=header,
                                     encoding=spec.get("encoding", "utf-8"), **options)
    elif suffix in {".sdf", ".mol"}:
        from .data_io import DataHandler
        smiles_col = spec.get("columns", {}).get("smiles", "smiles")
        frame = DataHandler._load_sdf(str(path), smiles_col=smiles_col)
    else:
        raise CliError("UNSUPPORTED_FORMAT", "Unsupported input format.", {"suffix": suffix})
    frame.columns = frame.columns.map(str)
    unit_column = spec.get("columns", {}).get("unit")
    if unit_column in frame:
        frame[unit_column] = frame[unit_column].astype("string")
    return frame


def inspect_file(spec, limit=5):
    path = Path(spec["path"])
    if not path.is_file():
        raise CliError("INPUT_NOT_FOUND", "Input file does not exist.", {"path": str(path)}, exit_code=3)
    sheets = _sheets(path) if path.suffix.lower() in {".xlsx", ".xls"} else []
    if len(sheets) > 1 and "sheet" not in spec:
        return envelope("inspect", summary={"path": str(path), "sheets": sheets, "requires_sheet": True})
    frame = read_frame(spec)
    return envelope("inspect", summary={
        "path": str(path), "sheets": sheets, "rows": len(frame),
        "columns": [{"name": name, "dtype": str(dtype)} for name, dtype in frame.dtypes.items()],
        "preview": json_safe(frame.head(limit).to_dict(orient="records")),
        "preview_limit": limit,
    })


def _check_outputs(config, protected_paths=()):
    outputs = {key: Path(value) for key, value in config["output"].items() if key.endswith("path")}
    protected = [Path(config["input"]["path"]), *map(Path, protected_paths)]
    seen = []
    for key, path in outputs.items():
        for other in protected + seen:
            same = path == other or (path.exists() and other.exists() and os.path.samefile(path, other))
            if same:
                raise CliError("OUTPUT_CONFLICT", "Output paths must be distinct from inputs, configuration and other outputs.", {"path": str(path)})
        if path.exists() and (not config["overwrite"] or not path.is_file()):
            raise CliError("OUTPUT_EXISTS", "Output exists; choose another path or use --overwrite.", {"path": str(path)}, exit_code=3)
        parent = path.parent
        while not parent.exists():
            parent = parent.parent
        if not parent.is_dir():
            raise CliError("INVALID_OUTPUT_DIRECTORY", "An output parent is not a directory.", {"path": str(parent)}, exit_code=3)
        seen.append(path)
    return outputs


def _unit_rules(params):
    from .unit_processor import UnitProcessor
    processor = UnitProcessor()
    for rule in params["conversion_rules"]:
        expression = rule["formula"].replace("^", "**")
        try:
            tree = ast.parse(expression, mode="eval")
        except SyntaxError as exc:
            raise CliError("INVALID_FORMULA", "Invalid unit conversion formula.", {"formula": rule["formula"]}) from exc
        allowed = (ast.Expression, ast.BinOp, ast.UnaryOp, ast.Add, ast.Sub, ast.Mult,
                   ast.Div, ast.Pow, ast.UAdd, ast.USub, ast.Call, ast.Name, ast.Load, ast.Constant)
        for node in ast.walk(tree):
            if not isinstance(node, allowed) or (isinstance(node, ast.Name) and node.id not in {"x", "mw", "log", "log10", "exp", "e"}):
                raise CliError("INVALID_FORMULA", "Unsupported expression in unit conversion formula.")
            if isinstance(node, ast.Constant) and (type(node.value) not in (int, float) or
                                                   (isinstance(node.value, float) and not math.isfinite(node.value))):
                raise CliError("INVALID_FORMULA", "Formula constants must be finite numbers.")
            if isinstance(node, ast.Call) and (not isinstance(node.func, ast.Name) or node.func.id not in {"log", "log10", "exp"} or len(node.args) != 1 or node.keywords):
                raise CliError("INVALID_FORMULA", "Only log(x), log10(x) and exp(x) calls are supported.")
        processor.add_rule(rule["from"], rule["to"], rule["formula"])
    return processor


def _preflight(config, frame):
    """Validate dependencies using input metadata, without running chemistry."""
    available = set(frame.columns)
    roles = dict(config["input"]["columns"])
    missing = sorted(set(roles.values()) - available)
    if missing:
        raise CliError("MISSING_COLUMNS", "Mapped columns were not found.", {"missing": missing, "available": list(frame.columns)})
    units = set(frame[roles["unit"]].dropna().astype(str)) - {""} if "unit" in roles else set()
    source_columns = list(frame.columns)
    aliases = {}

    def protect_columns(output_columns):
        renamed = source_column_renames(available, source_columns, output_columns)
        available.difference_update(renamed)
        available.update(renamed.values())
        source_columns[:] = [renamed.get(column, column) for column in source_columns]
        roles.update({role: renamed.get(column, column) for role, column in roles.items()})
        aliases.update({name: renamed.get(column, column) for name, column in aliases.items()})
        aliases.update(renamed)
        return renamed

    protect_columns(INITIAL_COLUMNS)
    available.update({SOURCE_ROWS, *INITIAL_COLUMNS})
    standardized = False
    preprocessed = False
    previous_units = None
    conversions = {}
    scales = {match.group(1) if (match := re.fullmatch(r"(-?log10)\(.*\)", unit)) else "linear"
              for unit in units} if "target" in roles else set()
    target_scale = next(iter(scales)) if len(scales) == 1 else "mixed" if scales else "linear"
    for number, step in enumerate(config["steps"], 1):
        params = step["params"]
        if step["op"] == "units":
            if target_scale != "linear":
                raise CliError("LOG_TRANSFORMED_TARGET", "Unit conversion requires a linear target; move the units step before the log transformation.",
                               {"step": number, "scale": target_scale})
            processor = _unit_rules(params)
            target = params["standard_unit"]
            current_roles = {role: roles[role] for role in ("target", "unit")}
            retry = (previous_units is not None
                     and target == previous_units["target"]
                     and current_roles == previous_units["output_roles"]
                     and set(previous_units["input_roles"].values()) <= available)
            input_roles = previous_units["input_roles"] if retry else current_roles
            input_units = previous_units["input_units"] if retry else units
            mw_column = aliases.get(params.get("molecular_weight_col"), params.get("molecular_weight_col"))
            mw_source = params.get("molecular_weight_source")
            requested_basis = f"column:{mw_column}" if mw_column else mw_source
            if requested_basis and units == {target}:
                cursor, visited = current_roles["target"], set()
                while cursor in conversions and cursor not in visited:
                    visited.add(cursor)
                    conversion = conversions[cursor]
                    if conversion["requires_mw"]:
                        if requested_basis != conversion["mw_basis"]:
                            input_roles = conversion["input_roles"]
                            input_units = conversion["input_units"]
                        break
                    cursor = conversion["input_roles"]["target"]
            missing_rules = [{"from": unit, "to": target} for unit in sorted(input_units)
                             if unit != target and not processor.get_rule(unit, target)]
            if missing_rules:
                raise CliError("MISSING_CONVERSION_RULE", "Provide conversion_rules for these unit pairs.", {"step": number, "pairs": missing_rules})
            if mw_column and mw_column not in available:
                raise CliError("MISSING_COLUMNS", "Molecular-weight column was not found.", {"missing": [mw_column], "step": number})
            uses_mw = any("mw" in (processor.get_rule(unit, target) or "") for unit in input_units if unit != target)
            if uses_mw and not mw_column and mw_source is None:
                raise CliError("MOLECULAR_WEIGHT_SOURCE_REQUIRED", "Choose molecular_weight_source ('original' or 'standardized') or molecular_weight_col for this conversion.", {"step": number})
            if mw_source == "standardized" and not preprocessed:
                raise CliError("PREPROCESSING_REQUIRED", "The standardized molecular-weight source requires a preceding preprocess step.", {"step": number})
            if uses_mw and not mw_column and "smiles" not in roles:
                raise CliError("MOLECULAR_WEIGHT_REQUIRED", "Map a SMILES column or specify molecular_weight_col.", {"step": number})
            while True:
                renamed = protect_columns({*(column + " (Standardized)" for column in input_roles.values()),
                                           "Unit Conversion Status", "Molecular Weight Used", "Molecular Weight Source"})
                if not renamed:
                    break
                input_roles = {role: renamed.get(column, column) for role, column in input_roles.items()}
            output_roles = {role: input_roles[role] + " (Standardized)" for role in ("target", "unit")}
            roles.update(output_roles)
            available.update([roles["target"], roles["unit"]])
            available.update({"Unit Conversion Status", "Molecular Weight Used", "Molecular Weight Source"})
            units = {target}
            previous_units = {"input_roles": dict(input_roles), "input_units": set(input_units),
                              "output_roles": output_roles, "target": target,
                              "requires_mw": uses_mw, "mw_basis": requested_basis}
            conversions[output_roles["target"]] = previous_units
            standardized = True
        elif step["op"] == "deduplicate":
            if target_scale == "mixed" and params["data_type"] != "smiles":
                raise CliError("LOG_TRANSFORMED_TARGET", "The target mixes linear and logarithmic scales; use values on one scale.",
                               {"step": number})
            missing = sorted({aliases.get(column, column) for column in params["condition_cols"]} - available)
            if missing:
                raise CliError("MISSING_COLUMNS", "Deduplication condition columns were not found.", {"missing": missing, "step": number})
            if len(units) > 1 and not standardized and params["data_type"] != "smiles":
                raise CliError("UNIT_STANDARDIZATION_REQUIRED", "Add a units step before deduplicating mixed units.", {"step": number, "units": sorted(units)})
            while True:
                output_columns = {"Deduplication Strategy", "Deduplication Record Count", "Deduplication Source Rows",
                                  "Deduplication Input Values", "Deduplication Distinct Value Count", "Deduplication Value Range"}
                if params["data_type"] != "smiles":
                    output_columns.add(roles["target"] + "_new")
                if not protect_columns(output_columns):
                    break
            if params["data_type"] != "smiles":
                output_scale = target_scale if params["data_type"] == "continuous" else "linear"
                transform = params["log_transform"]
                if params["data_type"] == "continuous" and transform != "None":
                    if target_scale not in {"linear", transform}:
                        raise CliError("LOG_TRANSFORMED_TARGET", "The target already uses a different logarithmic scale; use the original linear target.",
                                       {"step": number, "scale": target_scale, "requested_scale": transform})
                    output_scale = transform
                roles["target"] += "_new"
                available.add(roles["target"])
                available.add("Deduplication Input Values")
                if output_scale != target_scale and "unit" in roles:
                    output_unit = roles["unit"] + "_new"
                    protect_columns({output_unit})
                    roles["unit"] = output_unit
                    available.add(output_unit)
                    units = {f"{output_scale}({unit})" for unit in units}
                target_scale = output_scale
            available.update({"Deduplication Strategy", "Deduplication Record Count", "Deduplication Source Rows",
                              "Deduplication Distinct Value Count", "Deduplication Value Range"})
        elif step["op"] == "preprocess":
            protect_columns(PREPROCESS_COLUMNS)
            available.update(PREPROCESS_COLUMNS)
            preprocessed = True
        elif step["op"] in {"search", "filter-atoms", "inchi"}:
            from .cli_operations import preflight_operation
            preflight_operation(step["op"], params)
            generated = (f"Substructure_{params['query_pattern']}" if step["op"] == "search" else
                         "InChI" if step["op"] == "inchi" else None)
            structure_column = "Canonical SMILES" if preprocessed else roles.get("smiles")
            if generated is not None and generated == structure_column:
                raise CliError("RESULT_COLUMN_CONFLICT", "The generated result column would overwrite the active SMILES column.",
                               {"step": number, "column": generated})
            if generated is not None:
                protect_columns({generated})
            if step["op"] == "search":
                available.add(f"Substructure_{params['query_pattern']}")
            elif step["op"] == "inchi":
                available.add("InChI")
        elif step["op"] == "enrich":
            from .cli_network import preflight_enrich
            preflight_enrich(params, roles)
            available.update(f"{prop}_from_web" for prop in params["request"])
            available.update({"Query_Status", "Query_Missing_Fields", "Query_Errors", "Data_Source", "Query_Method"})
            for role in ("smiles", "cas", "name"):
                if role in params["request"] and not (role == "smiles" and preprocessed):
                    roles[role] = f"{role}_from_web"
    missing = sorted(set(config["output"].get("columns", [])) - available)
    if missing:
        raise CliError("MISSING_COLUMNS", "Output columns will not be produced by this pipeline.", {"missing": missing})
    if Path(config["output"]["path"]).suffix.lower() in {".sdf", ".smi"}:
        structure_column = "Canonical SMILES" if preprocessed else roles.get("smiles")
        selected = set(config["output"].get("columns") or available)
        if not structure_column or structure_column not in selected:
            raise CliError("STRUCTURE_COLUMN_REQUIRED", "Structure export requires the active SMILES column in output.columns.",
                           {"column": structure_column})


def _replace_nonfinite_targets(pipeline):
    """Treat infinities as missing values before numeric algorithms run."""
    import pandas as pd
    if pipeline.target_col:
        numeric = pd.to_numeric(pipeline.df[pipeline.target_col], errors="coerce")
        mask = numeric.isin([float("inf"), float("-inf")])
        if mask.any():
            pipeline.df.loc[mask, pipeline.target_col] = pd.NA


def _audit_frame(frame, step, operation, reasons, severity="error"):
    result = frame.copy()
    result["CLI Step"] = step
    result["CLI Operation"] = operation
    result["CLI Reason"] = reasons
    result["CLI Severity"] = severity
    return result


def _write_frame(frame, path, smiles_col, id_col=None, source_smiles_col=None):
    from .data_io import DataHandler
    if path.suffix.lower() in {".sdf", ".smi"}:
        if not smiles_col or smiles_col not in frame:
            raise CliError("STRUCTURE_COLUMN_REQUIRED", "Structure export requires a mapped SMILES column.")
        if path.suffix.lower() == ".smi":
            from rdkit import Chem
            invalid = frame[smiles_col].map(lambda value: not isinstance(value, str) or Chem.MolFromSmiles(value) is None)
            if invalid.any():
                raise CliError("INVALID_STRUCTURE_EXPORT", "Use CSV/XLSX/TXT or SDF to retain invalid structures, or deduplicate valid results before SMI export.", {"invalid_rows": int(invalid.sum())})
    DataHandler.save_data(frame.copy(), str(path), frame.columns.tolist(), smiles_col, id_col,
                          interactive=False, source_smiles_col=source_smiles_col)


def _publish(config, outputs, pipeline, excluded, report):
    """Stage every artifact first; each destination is committed atomically."""
    staged = {}
    committed = {}
    previous_report = None
    try:
        for key, destination in outputs.items():
            destination.parent.mkdir(parents=True, exist_ok=True)
            descriptor, temporary = tempfile.mkstemp(prefix=".diptox-", suffix=destination.suffix.lower(), dir=destination.parent)
            os.close(descriptor)
            path = Path(temporary)
            staged[key] = path
            if key == "report_path":
                path.write_text(json.dumps(json_safe(report), ensure_ascii=False, allow_nan=False, indent=2) + "\n", encoding="utf-8")
            elif key == "excluded_path":
                _write_frame(excluded, path, pipeline.smiles_col, pipeline.id_col, pipeline.smiles_col)
            else:
                columns = config["output"].get("columns")
                frame = pipeline.df
                if columns:
                    missing = sorted(set(columns) - set(frame.columns))
                    if missing:
                        raise CliError("MISSING_COLUMNS", "Output columns were not found.", {"missing": missing})
                    frame = frame[columns]
                smiles = "Canonical SMILES" if pipeline._preprocess_key else pipeline.smiles_col
                _write_frame(frame, path, smiles, pipeline.id_col, pipeline.smiles_col)
        # The completion report is always published last.
        report_path = outputs.get("report_path")
        if config["overwrite"] and report_path and report_path.exists():
            descriptor, backup = tempfile.mkstemp(prefix=".diptox-previous-", suffix=".json", dir=report_path.parent)
            os.close(descriptor)
            previous_report = Path(backup)
            # A previous success report must not describe partially new data.
            os.replace(report_path, previous_report)
        order = sorted(outputs, key=lambda key: key == "report_path")
        for key in order:
            destination = outputs[key]
            if config["overwrite"]:
                os.replace(staged[key], destination)
            else:
                # Atomic no-clobber creation, including a concurrent writer.
                os.link(staged[key], destination)
                staged[key].unlink()
            committed[key] = str(destination)
    except BaseException as exc:
        if isinstance(exc, CliError):
            raise
        if isinstance(exc, (OSError, KeyboardInterrupt)):
            raise CliError("CANCELLED" if isinstance(exc, KeyboardInterrupt) else "WRITE_FAILED",
                           "Artifact publication did not complete.",
                           {"reason": str(exc), "published_artifacts": committed},
                           exit_code=130 if isinstance(exc, KeyboardInterrupt) else 3) from exc
        raise
    finally:
        if previous_report is not None:
            if not committed and not outputs["report_path"].exists():
                os.replace(previous_report, outputs["report_path"])
            else:
                previous_report.unlink(missing_ok=True)
        for temporary in staged.values():
            temporary.unlink(missing_ok=True)


def execute(config, command="run", dry_run=False, protected_paths=()):
    import pandas as pd
    from .core import DiptoxPipeline
    started = time.perf_counter()
    outputs = _check_outputs(config, protected_paths)
    frame = read_frame(config["input"])
    pipeline = DiptoxPipeline(interactive=False)
    from .cli_operations import apply_rules
    effective_rules = apply_rules(pipeline.chem_processor, config.get("rules", {}))
    preflight_frame = frame
    if 'Structure Reliability' in frame:
        preflight_frame = frame.loc[~frame['Structure Reliability'].eq('Unreliable').fillna(False)]
    _preflight(config, preflight_frame)
    report = envelope(command, effective_config=config, processing_rules={
        "sha256": hashlib.sha256(json.dumps(effective_rules, sort_keys=True).encode("utf-8")).hexdigest(),
        "counts": {key: len(value) for key, value in effective_rules.items()},
    })
    if dry_run:
        report.update(dry_run=True, summary={"input_rows": len(frame), "steps": len(config["steps"]), "validation": "configuration and input metadata"})
        return report
    if SOURCE_ROWS not in frame:
        frame[SOURCE_ROWS] = [str(index) for index in range(len(frame))]
    else:
        frame[SOURCE_ROWS] = frame[SOURCE_ROWS].astype("string")
        if not frame[SOURCE_ROWS].str.fullmatch(r"\d+(;\d+)*").fillna(False).all():
            raise CliError("INVALID_PROVENANCE", f"Reserved column '{SOURCE_ROWS}' must contain semicolon-separated row numbers.")
    pipeline.load_data(frame, **{role + "_col": column for role, column in config["input"]["columns"].items()})
    audits = []
    steps = []
    warnings = []
    if not pipeline.excluded_df.empty:
        rejected = pipeline.excluded_df
        audits.append(_audit_frame(rejected, 0, 'load', rejected['Exclusion Reason']))
        warnings.append({'code': 'UNRELIABLE_STRUCTURE', 'step': 0, 'operation': 'load', 'rows': len(rejected)})
        steps.append({'step': 0, 'operation': 'load', 'input_rows': len(frame),
                      'output_rows': len(pipeline.df), 'rejected_rows': len(rejected),
                      'filtered_rows': 0, 'consolidated_rows': 0})
    for number, step in enumerate(config["steps"], 1):
        operation = step["op"]
        params = dict(step["params"])
        before = pipeline.df.copy()
        merged = 0
        operation_stats = {}
        if operation == "preprocess":
            pipeline.preprocess(**params)
            rejected = pipeline.df.loc[~pipeline.df["Is Valid"].fillna(False)]
            audit = _audit_frame(rejected, number, operation, rejected["Processing Log"])
        elif operation == "units":
            _replace_nonfinite_targets(pipeline)
            params["conversion_rules"] = {(rule["from"], rule["to"]): rule["formula"] for rule in params["conversion_rules"]}
            pipeline.standardize_units(**params)
            nonfinite = pd.to_numeric(pipeline.df[pipeline.target_col], errors="coerce").isin([float("inf"), float("-inf")])
            if nonfinite.any():
                pipeline.df.loc[nonfinite, pipeline.target_col] = pd.NA
                pipeline.df.loc[nonfinite, pipeline.unit_col] = pd.NA
                pipeline.df.loc[nonfinite, "Unit Conversion Status"] = "Non-finite conversion result"
            rejected = pipeline.df.loc[~pipeline.df["Unit Conversion Status"].isin(["Converted", "No conversion needed"])]
            audit = _audit_frame(rejected, number, operation, rejected["Unit Conversion Status"])
        elif operation == "deduplicate":
            # A structure-only operation must not acquire target/unit semantics.
            original_target, original_unit = pipeline.target_col, pipeline.unit_col
            if params["data_type"] == "smiles":
                pipeline.target_col = pipeline.unit_col = None
            elif params["data_type"] == "continuous":
                _replace_nonfinite_targets(pipeline)
            pipeline.config_deduplicator(**params)
            pipeline.dataset_deduplicate()
            if params["data_type"] == "continuous":
                aggregate = pd.to_numeric(pipeline.df[pipeline.target_col], errors="coerce")
                if aggregate.isin([float("inf"), float("-inf")]).any():
                    raise CliError("NONFINITE_AGGREGATE", "Deduplication produced a non-finite value; check the input scale.", {"step": number})
            if params["data_type"] == "smiles":
                pipeline.target_col, pipeline.unit_col = original_target, original_unit
            reasons = pipeline.deduplicator.exclusion_reasons
            rejected = before.loc[list(reasons)]
            audit = _audit_frame(rejected, number, operation, [reasons[index] for index in rejected.index])
            merged = max(0, len(before) - len(pipeline.df) - len(rejected))
            if not pipeline.df.empty:
                provenance = before[SOURCE_ROWS].to_dict()
                pipeline.df[SOURCE_ROWS] = [
                    ';'.join(dict.fromkeys(
                        source for member in pipeline.deduplicator.source_indices[index]
                        for source in str(provenance[member]).split(';')
                    )) for index in pipeline.df.index
                ]
        elif operation in {"search", "filter-atoms", "inchi"}:
            from .cli_operations import execute_operation
            audit, operation_stats = execute_operation(pipeline, operation, params)
            audit = _audit_frame(audit, number, operation, audit["CLI Reason"], audit["CLI Severity"])
        elif operation == "enrich":
            from .cli_network import enrich_pipeline
            operation_stats = enrich_pipeline(pipeline, params)
            incomplete = pipeline.df.loc[pipeline.df["Query_Status"] != "complete"]
            reasons = "Missing fields: " + incomplete["Query_Missing_Fields"].astype(str) + " | " + incomplete["Query_Errors"].astype(str)
            audit = _audit_frame(incomplete, number, operation, reasons)
        else:
            raise CliError("INVALID_OPERATION", f"Unknown operation: {operation}")
        rejected_count = int(audit["CLI Severity"].ne("filtered").sum())
        filtered_count = int(audit["CLI Severity"].eq("filtered").sum())
        if not audit.empty:
            audits.append(audit)
        if rejected_count:
            warning_code = {"units": "UNIT_CONVERSION_INCOMPLETE", "enrich": "ENRICHMENT_INCOMPLETE"}.get(operation, "RECORDS_REJECTED")
            warnings.append({"code": warning_code, "step": number, "operation": operation, "rows": rejected_count})
        steps.append({"step": number, "operation": operation, "input_rows": len(before),
                      "output_rows": len(pipeline.df), **operation_stats, "rejected_rows": rejected_count,
                      "filtered_rows": filtered_count, "consolidated_rows": merged})
    excluded = pd.concat(audits, ignore_index=True) if audits else pd.DataFrame(columns=list(dict.fromkeys([*frame.columns, *AUDIT_COLUMNS])))
    if config["strict"] and warnings:
        raise CliError("RECORD_FAILURES", "Strict mode rejected record-level processing failures; no artifacts were written.", {"steps": steps, "warnings": warnings}, exit_code=5)
    report.update(
        summary={"input_rows": len(frame), "output_rows": len(pipeline.df), "exclusion_events": len(excluded),
                 "steps": steps, "columns": list(pipeline.df.columns),
                 "column_roles": {"smiles": "Canonical SMILES" if pipeline._preprocess_key else pipeline.smiles_col,
                                  "target": pipeline.target_col, "unit": pipeline.unit_col,
                                  "cas": pipeline.cas_col, "name": pipeline.name_col,
                                  "id": pipeline.id_col, "inchikey": pipeline.inchikey_col}},
        warnings=warnings, artifacts={key: str(path) for key, path in outputs.items()},
        history=json_safe(pipeline.get_history().to_dict(orient="records")),
        completed_at=datetime.now(timezone.utc).isoformat(),
        elapsed_seconds=round(time.perf_counter() - started, 6),
        runtime={"python": platform.python_version(), "pandas": pd.__version__},
    )
    from rdkit import __version__ as rdkit_version
    report["runtime"]["rdkit"] = rdkit_version
    _publish(config, outputs, pipeline, excluded, report)
    return report
