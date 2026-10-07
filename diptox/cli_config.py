"""Dependency-free configuration and discovery contract for the DiPTox CLI."""

from copy import deepcopy
import math
from pathlib import Path


class CliError(Exception):
    """An expected CLI failure with a stable machine-readable identity."""

    def __init__(self, code: str, message: str, details=None,
                 retryable: bool = False, exit_code: int = 2):
        super().__init__(message)
        self.code = code
        self.message = message
        self.details = {} if details is None else details
        self.retryable = retryable
        self.exit_code = exit_code


INPUT_FORMATS = (".csv", ".xlsx", ".xls", ".txt", ".smi", ".sdf", ".mol")
OUTPUT_FORMATS = (".csv", ".xlsx", ".txt", ".smi", ".sdf")
COLUMN_ROLES = ("smiles", "target", "unit", "cas", "name", "inchikey", "id")
NETWORK_SOURCES = ("pubchem", "chemspider", "comptox", "cactus", "chembl", "cas")
NETWORK_API_KEY_ENV = {
    "chemspider": "DIPTOX_CHEMSPIDER_API_KEY",
    "comptox": "DIPTOX_COMPTOX_API_KEY",
    "cas": "DIPTOX_CAS_API_KEY",
}
_STRING = {"type": "string", "minLength": 1}
_STRINGS = {"type": "array", "items": _STRING, "default": []}
_BOOL_DEFAULTS = {
    "remove_salts": True,
    "remove_solvents": True,
    "remove_inorganic": True,
    "neutralize": True,
    "reject_non_neutral": False,
    "remove_stereo": False,
    "remove_isotopes": True,
    "remove_hs": True,
    "add_hs": False,
    "sanitize": True,
    "reject_radical_species": True,
}

PARAM_SPECS = {
    "preprocess": {
        "mixture_mode": {"type": "string", "enum": ["keep", "reject", "largest"],
                         "default": "reject", "description": "Default: reject. largest rejects tied largest fragments as Ambiguous parent."},
        "element_policy": {"type": "string",
                           "enum": ["allow_all", "reject_metals", "allowed_atoms"],
                           "default": "allow_all"},
        **{name: {"type": "boolean", "default": default}
           for name, default in _BOOL_DEFAULTS.items()},
        "hac_threshold": {"type": "integer", "minimum": 0, "default": 3},
        "n_jobs": {"type": "integer", "default": 1,
                   "anyOf": [{"const": -1}, {"minimum": 1}]},
        "chunksize": {"type": "integer", "minimum": 1, "default": 100},
    },
    "units": {
        "standard_unit": {"type": "string", "minLength": 1},
        "molecular_weight_source": {"type": ["string", "null"],
                                    "enum": [None, "original", "standardized"],
                                    "default": None,
                                    "description": "Required for MW-dependent conversion unless molecular_weight_col is provided. Choose original or standardized; standardized requires preprocessing."},
        "molecular_weight_col": {"type": ["string", "null"], "minLength": 1,
                                 "default": None},
        "conversion_rules": {
            "type": "array", "default": [],
            "items": {"type": "object", "additionalProperties": False,
                      "required": ["from", "to", "formula"],
                      "properties": {name: _STRING for name in ("from", "to", "formula")}},
        },
    },
    "deduplicate": {
        "data_type": {"type": "string", "enum": ["smiles", "discrete", "continuous"],
                      "default": "continuous"},
        "method": {"type": "string", "enum": ["auto", "IQR", "3sigma",
                                               "vote", "priority"],
                   "default": "auto", "description": "Continuous outlier filtering; discrete selection defaults to vote."},
        "aggregation": {"type": "string", "enum": ["mean", "max", "min"],
                        "default": "mean", "description": "Continuous aggregation after target transformation and outlier filtering."},
        "condition_cols": _STRINGS,
        "priority": {"type": "array", "default": [],
                     "items": {"type": ["string", "number", "boolean", "null"]}},
        "p_threshold": {"type": "number", "exclusiveMinimum": 0,
                        "exclusiveMaximum": 1, "default": 0.05},
        "log_transform": {"type": "string", "enum": ["None", "-log10", "log10"],
                          "default": "None"},
        "dropna_conditions": {"type": "boolean", "default": False},
    },
    "search": {
        "query_pattern": _STRING,
        "query_type": {"type": "string", "enum": ["smarts", "smiles"], "default": "smarts"},
        "mode": {"type": "string", "enum": ["annotate", "matches", "nonmatches"],
                 "default": "annotate"},
    },
    "filter-atoms": {
        name: {"type": ["integer", "null"], "minimum": 0, "default": None}
        for name in ("min_heavy_atoms", "max_heavy_atoms", "min_total_atoms", "max_total_atoms")
    },
    "filter-values": {
        "column": _STRING,
        "mode": {"type": "string", "enum": ["keep", "remove"], "default": "keep"},
        "values": {"type": "array", "default": [],
                   "items": {"type": ["string", "number", "boolean", "null"]},
                   "description": "JSON array of values; [] keeps all rows, null selects missing cells."},
    },
    "inchi": {},
    "enrich": {
        "sources": {"type": "array", "minItems": 1, "uniqueItems": True,
                    "items": {"type": "string", "enum": list(NETWORK_SOURCES)},
                    "default": ["pubchem"]},
        "send": {"type": "array", "minItems": 1, "uniqueItems": True,
                 "items": {"type": "string", "enum": ["smiles", "cas", "name"]},
                 "default": ["smiles"]},
        "request": {"type": "array", "minItems": 1, "uniqueItems": True,
                    "items": {"type": "string", "enum": ["smiles", "cas", "iupac", "mw", "name"]}},
        "max_workers": {"type": "integer", "minimum": 1, "maximum": 32, "default": 4},
        "timeout": {"type": "number", "exclusiveMinimum": 0, "default": 15},
        "deadline": {"type": "number", "exclusiveMinimum": 0, "default": 120},
        "retries": {"type": "integer", "minimum": 0, "maximum": 10, "default": 2,
                    "description": "Additional attempts after the initial request."},
        "retry_delay": {"type": "number", "minimum": 0, "default": 1},
        "interval": {"type": "number", "minimum": 0, "default": 0.3},
    },
}

PARAM_SPECS['transform-column'] = {
    **PARAM_SPECS['units'],
    'value_col': _STRING,
    'unit_col': _STRING,
    'standard_unit': {'type': ['string', 'null'], 'minLength': 1, 'default': None},
    'log_transform': PARAM_SPECS['deduplicate']['log_transform'],
}

PARAM_SPECS['merge-values'] = {
    'groups': {'type': ['array', 'null'], 'default': None,
               'items': {'type': 'object', 'additionalProperties': False,
                         'required': ['values', 'replacement'],
                         'properties': {'values': {'type': 'array', 'minItems': 1, 'items': {'type': ['string', 'number', 'boolean', 'null']}},
                                        'replacement': {'type': 'string', 'minLength': 1}}}},
    'column': _STRING,
    'values': PARAM_SPECS['filter-values']['values'],
    'replacement': {'type': 'string', 'minLength': 1, 'default': 'other'},
    'output_column': {'type': ['string', 'null'], 'minLength': 1, 'default': None},
}

REQUIRED_PARAMS = {
    op: {"units": ("standard_unit",), "search": ("query_pattern",),
         "enrich": ("request",), "filter-values": ("column",),
         "transform-column": ("value_col", "unit_col"), "merge-values": ("column",)}.get(op, ())
    for op in PARAM_SPECS
}


def _object(properties, required=()):
    return {"type": "object", "properties": properties,
            "required": list(required), "additionalProperties": False}


_RULE_CHANGES = {**_object({
    "add": {**_STRINGS, "uniqueItems": True},
    "remove": {**_STRINGS, "uniqueItems": True},
}), "default": {}}
RULES_SCHEMA = {
    **_object({
        "atoms": _RULE_CHANGES,
        "salts": _RULE_CHANGES,
        "solvents": _RULE_CHANGES,
        "neutralization": {**_object({
            "add": {"type": "array", "default": [], "uniqueItems": True,
                    "items": _object({"reactant": _STRING, "product": _STRING},
                                     ("reactant", "product"))},
            "remove": {**_STRINGS, "uniqueItems": True},
        }), "default": {}},
    }),
    "default": {},
}


def _config_schema():
    step_schemas = []
    for op, params in PARAM_SPECS.items():
        required = REQUIRED_PARAMS[op]
        params_schema = _object(params, required)
        if not required:
            params_schema["default"] = {}
        step_schemas.append(_object({"op": {"const": op}, "params": params_schema},
                                    ("op", "params") if required else ("op",)))
    return _object({
        "schema_version": {"type": "string", "const": "1"},
        "input": _object({
            "path": _STRING,
            "columns": {**_object({role: _STRING for role in COLUMN_ROLES}), "default": {}},
            "sheet": {"type": ["string", "integer"], "minLength": 1, "minimum": 0},
            "header": {"type": "boolean", "default": True},
            "delimiter": {"type": "string", "minLength": 1},
            "encoding": {"type": "string", "minLength": 1, "default": "utf-8"},
        }, ("path",)),
        "steps": {"type": "array", "minItems": 1, "items": {"oneOf": step_schemas}},
        "output": _object({
            "path": _STRING,
            "excluded_path": _STRING,
            "report_path": _STRING,
            "columns": {"type": "array", "minItems": 1, "items": _STRING},
        }, ("path",)),
        "overwrite": {"type": "boolean", "default": False},
        "strict": {"type": "boolean", "default": False},
        "rules": RULES_SCHEMA,
    }, ("schema_version", "input", "steps", "output"))


def _invalid(field, message, **details):
    raise CliError("INVALID_CONFIG", f"{field}: {message}", {"field": field, **details})


def _validate(value, spec, field):
    """Validate the small JSON-schema vocabulary used by this module."""
    expected = spec.get("type")
    expected = [expected] if isinstance(expected, str) else expected
    types = {
        "object": isinstance(value, dict), "array": isinstance(value, list),
        "string": isinstance(value, str), "boolean": type(value) is bool,
        "integer": type(value) is int,
        "number": type(value) is int or (type(value) is float and math.isfinite(value)),
        "null": value is None,
    }
    if expected and not any(types.get(kind, False) for kind in expected):
        _invalid(field, "has an invalid type", expected=expected)
    if type(value) is float and not math.isfinite(value):
        _invalid(field, "must be a finite number")
    if "const" in spec and value != spec["const"]:
        _invalid(field, f"must equal {spec['const']!r}")
    if "enum" in spec and value not in spec["enum"]:
        _invalid(field, "is not an allowed value", allowed=spec["enum"])
    if isinstance(value, str) and spec.get("minLength") and not value.strip():
        _invalid(field, "must not be empty or whitespace")
    if type(value) in (int, float):
        checks = (("minimum", lambda a, b: a >= b), ("maximum", lambda a, b: a <= b),
                  ("exclusiveMinimum", lambda a, b: a > b),
                  ("exclusiveMaximum", lambda a, b: a < b))
        for bound, compare in checks:
            if bound in spec and not compare(value, spec[bound]):
                _invalid(field, f"does not satisfy {bound}={spec[bound]}")
    if "anyOf" in spec:
        for option in spec["anyOf"]:
            try:
                _validate(value, option, field)
                break
            except CliError:
                pass
        else:
            _invalid(field, "does not match any allowed value range")
    if isinstance(value, dict) and "properties" in spec:
        properties = spec["properties"]
        unknown = [key for key in value if key not in properties]
        if unknown:
            _invalid(field, "contains unknown fields", unknown_fields=unknown)
        missing = [key for key in spec.get("required", []) if key not in value]
        if missing:
            _invalid(field, "is missing required fields", missing_fields=missing)
        result = {}
        for name, child_spec in properties.items():
            if name in value:
                result[name] = _validate(value[name], child_spec, f"{field}.{name}")
            elif "default" in child_spec:
                result[name] = _validate(deepcopy(child_spec["default"]), child_spec, f"{field}.{name}")
        return result
    if isinstance(value, list):
        if len(value) < spec.get("minItems", 0):
            _invalid(field, "must not be empty")
        if spec.get("uniqueItems") and any(item in value[:index] for index, item in enumerate(value)):
            _invalid(field, "must not contain duplicate values")
        return [_validate(item, spec.get("items", {}), f"{field}[{index}]")
                for index, item in enumerate(value)]
    return value


def normalize_rules(raw: dict) -> dict:
    """Validate custom rules structurally without importing chemistry dependencies."""
    return _validate(raw, RULES_SCHEMA, "rules")


def _resolve_path(value, base_dir, field, formats=None):
    try:
        path = Path(value).expanduser()
        resolved = (path if path.is_absolute() else base_dir / path).resolve()
    except (OSError, ValueError, RuntimeError) as error:
        _invalid(field, f"invalid path: {error}")
    if formats and resolved.suffix.lower() not in formats:
        _invalid(field, "unsupported file extension", allowed_extensions=list(formats))
    return str(resolved)


def normalize_config(raw: dict, base_dir: Path) -> dict:
    """Validate JSON configuration, add defaults, and resolve its file paths.

    Relative paths are based on the configuration file's directory. This function
    never opens input files or creates directories; execution performs those checks.
    """
    schema = _config_schema()
    # Each operation is validated below so errors identify its exact field.
    schema["properties"]["steps"]["items"] = _object({"op": {"type": "string",
                                                               "enum": list(PARAM_SPECS)},
                                                       "params": {"type": "object", "default": {}}},
                                                      ("op",))
    config = _validate(raw, schema, "config")
    available_roles = set(config["input"]["columns"])
    for index, step in enumerate(config["steps"]):
        op = step["op"]
        field = f"config.steps[{index}].params"
        params = step["params"]
        method_given = "method" in params
        params = _validate(params, _object(PARAM_SPECS[op], REQUIRED_PARAMS[op]), field)
        if op == "deduplicate":
            data_type = params["data_type"]
            if data_type == "discrete" and not method_given:
                params["method"] = "vote"
            methods = {"continuous": ["auto", "IQR", "3sigma"],
                       "discrete": ["vote", "priority"], "smiles": ["auto"]}
            if params["method"] not in methods[data_type]:
                _invalid(f"{field}.method", f"is incompatible with data_type={data_type!r}",
                         allowed=methods[data_type])
            if data_type != "continuous" and params["aggregation"] != "mean":
                _invalid(f"{field}.aggregation", "requires continuous data")
            if params["method"] == "priority" and not params["priority"]:
                _invalid(f"{field}.priority", "must contain ordered values when method is priority")
            if params["priority"] and params["method"] != "priority":
                _invalid(f"{field}.priority", "is only supported by the priority method")
            if params["log_transform"] != "None" and data_type != "continuous":
                _invalid(f"{field}.log_transform", "requires continuous data")
        if op in {"units", "transform-column"}:
            seen = set()
            for rule in params["conversion_rules"]:
                key = (rule["from"], rule["to"])
                if key in seen:
                    _invalid(f"{field}.conversion_rules", "contains a duplicate from/to pair",
                             from_unit=key[0], to_unit=key[1])
                seen.add(key)
        if op == "filter-atoms":
            if all(value is None for value in params.values()):
                _invalid(field, "requires at least one atom-count bound")
            for kind in ("heavy", "total"):
                lower, upper = params[f"min_{kind}_atoms"], params[f"max_{kind}_atoms"]
                if lower is not None and upper is not None and lower > upper:
                    _invalid(field, f"min_{kind}_atoms must not exceed max_{kind}_atoms")
        required_roles = ([] if op in {"filter-values", "transform-column", "merge-values"} else ["target", "unit"] if op == "units" else
                          params["send"] if op == "enrich" else ["smiles"])
        if op == "deduplicate" and params["data_type"] != "smiles":
            required_roles.append("target")
        missing_roles = [role for role in required_roles if role not in available_roles]
        if missing_roles:
            _invalid("config.input.columns", f"is missing mappings required by step {index} ({op})",
                     missing_fields=missing_roles)
        step["params"] = params
        if op == "enrich":
            available_roles.update(set(params["request"]) & set(COLUMN_ROLES))

    base_dir = Path(base_dir)
    config["input"]["path"] = _resolve_path(config["input"]["path"], base_dir,
                                            "config.input.path", INPUT_FORMATS)
    output = config["output"]
    output["path"] = _resolve_path(output["path"], base_dir, "config.output.path", OUTPUT_FORMATS)
    if "excluded_path" in output:
        output["excluded_path"] = _resolve_path(output["excluded_path"], base_dir,
                                                "config.output.excluded_path", (".csv", ".xlsx", ".txt"))
    if "report_path" in output:
        output["report_path"] = _resolve_path(output["report_path"], base_dir,
                                              "config.output.report_path", (".json",))
    return config


def schema_document() -> dict:
    """Return a fresh, serializable discovery document without heavy imports."""
    return deepcopy({
        "schema_version": "1",
        "commands": ["schema", "inspect", "column-values", "run", *PARAM_SPECS, "rules"],
        "config_schema": {"$schema": "https://json-schema.org/draft/2020-12/schema",
                          **_config_schema()},
        "step_parameters": PARAM_SPECS,
        "required_step_parameters": REQUIRED_PARAMS,
        "rules_schema": RULES_SCHEMA,
        "network": {
            "sources": list(NETWORK_SOURCES),
            "api_key_environment": NETWORK_API_KEY_ENV,
            "credentials": "API keys are read from environment variables, never from pipeline configuration.",
        },
        "input_formats": list(INPUT_FORMATS),
        "output_formats": list(OUTPUT_FORMATS),
        "column_roles": list(COLUMN_ROLES),
        "path_resolution": "Relative to the configuration file directory; command-line paths use the working directory.",
        "continuous_aggregations": ["mean", "max", "min"],
        "deduplication_methods": {"continuous": ["auto", "IQR", "3sigma"],
                                  "discrete": ["vote", "priority"], "smiles": ["auto"]},
    })
