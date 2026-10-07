"""Non-interactive, JSON-first command line interface for DiPTox."""

from __future__ import annotations

import argparse
from contextlib import redirect_stdout
import json
import logging
from pathlib import Path
import sys

from . import __version__
from .cli_config import COLUMN_ROLES, INPUT_FORMATS, PARAM_SPECS, CliError, normalize_config, normalize_rules, schema_document


try:
    from argparse import BooleanOptionalAction
except ImportError:  # Python 3.8: preserve --flag / --no-flag and omitted values.
    class BooleanOptionalAction(argparse.Action):
        def __init__(self, option_strings, dest, default=None, **kwargs):
            paired_options = []
            for option in option_strings:
                paired_options.append(option)
                if option.startswith("--"):
                    paired_options.append("--no-" + option[2:])
            super().__init__(paired_options, dest, nargs=0, default=default, **kwargs)

        def __call__(self, parser, namespace, values, option_string=None):
            setattr(namespace, self.dest, not option_string.startswith("--no-"))

        def format_usage(self):
            return " | ".join(self.option_strings)


class Parser(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("allow_abbrev", False)
        super().__init__(*args, **kwargs)

    def error(self, message):
        raise CliError("INVALID_ARGUMENTS", message)


def _output_options(parser):
    parser.add_argument("--format", choices=["json", "text"], default=argparse.SUPPRESS,
                        help="Response format (default: json); does not change exported file format.")
    parser.add_argument("--verbose", action="store_true", default=argparse.SUPPRESS, help="Write informational logs to stderr.")
    parser.add_argument("--debug", action="store_true", default=argparse.SUPPRESS, help="Include a traceback on stderr for failures.")


def _input_options(parser):
    parser.add_argument("--input", required=True, help="Input file path.")
    sheets = parser.add_mutually_exclusive_group()
    sheets.add_argument("--sheet", default=argparse.SUPPRESS, help="Excel worksheet name (including numeric names).")
    sheets.add_argument("--sheet-index", type=int, default=argparse.SUPPRESS, help="Zero-based Excel worksheet index.")
    parser.add_argument("--header", action=BooleanOptionalAction, default=argparse.SUPPRESS,
                        help="Use --no-header for headerless files; column names then are 0, 1, ...")
    parser.add_argument("--delimiter", default=argparse.SUPPRESS, help="CSV/TXT delimiter.")
    parser.add_argument("--encoding", default=argparse.SUPPRESS, help="Text input encoding (default: utf-8).")
    for role in COLUMN_ROLES:
        parser.add_argument("--" + role + "-col", default=argparse.SUPPRESS, help=f"Column mapped to {role}.")


def build_parser():
    parser = Parser(prog="diptox", description="Non-interactive molecular data processing. Results are JSON; logs use stderr.")
    parser.add_argument("--version", action="version", version=f"diptox {__version__}")
    _output_options(parser)
    parser.set_defaults(format="json", verbose=False, debug=False)
    commands = parser.add_subparsers(dest="command", required=True, parser_class=Parser)
    schema = commands.add_parser("schema", help="Discover commands, configuration schema and defaults without loading chemistry libraries.")
    _output_options(schema)
    rules = commands.add_parser("rules", help="Show effective chemical rules, optionally applying a JSON rules file.")
    _output_options(rules)
    rules.add_argument("--rules-file", help="JSON rule changes to validate and apply for this invocation only.")
    inspect = commands.add_parser("inspect", help="List sheets, columns and a bounded preview (loads the selected table).")
    _output_options(inspect)
    _input_options(inspect)
    inspect.add_argument("--limit", type=int, default=5, help="Preview rows, 0 to 1000 (default: 5).")
    values = commands.add_parser("column-values", help="List every distinct value and row count in a column.")
    _output_options(values)
    _input_options(values)
    values.add_argument("--column", required=True)
    run = commands.add_parser("run", help="Execute ordered pipeline steps from JSON.")
    _output_options(run)
    run.add_argument("--config", required=True, help="JSON file, or - to read configuration from stdin.")
    run.add_argument("--rules-file", help="JSON rule changes; replaces the configuration's rules block.")
    for name in ("dry-run", "overwrite", "strict"):
        run.add_argument("--" + name, action="store_true", default=argparse.SUPPRESS)
    for operation, specs in PARAM_SPECS.items():
        command = commands.add_parser(operation, help=f"Run the {operation} step on a file.")
        _output_options(command)
        _input_options(command)
        command.add_argument("--output", required=True, help="Output file; suffix selects data format.")
        command.add_argument("--excluded", default=argparse.SUPPRESS, help="CSV/XLSX/TXT file with record-level audit events.")
        command.add_argument("--report", default=argparse.SUPPRESS, help="JSON completion report path.")
        command.add_argument("--columns", nargs="+", default=argparse.SUPPRESS, help="Export only these result columns.")
        command.add_argument("--rules-file", help="JSON chemical rule changes for this invocation only.")
        for name in ("dry-run", "overwrite", "strict"):
            command.add_argument("--" + name, action="store_true", default=argparse.SUPPRESS)
        for name, spec in specs.items():
            if operation == 'transform-column' and name == 'unit_col':
                command.add_argument('--column-unit', default=argparse.SUPPRESS,
                                     help='Unit column paired with --value-col (independent of the primary target).')
                continue
            options = {"default": argparse.SUPPRESS}
            kind = spec["type"]
            if isinstance(kind, list):
                kind = next(item for item in kind if item != "null")
            if kind == "boolean":
                options["action"] = BooleanOptionalAction
            elif kind == "integer":
                options["type"] = int
            elif kind == "number":
                options["type"] = float
            elif kind == "array" and name not in {"conversion_rules", "values", "groups"}:
                options["nargs"] = "+"
                if "enum" in spec.get("items", {}):
                    options["choices"] = spec["items"]["enum"]
            if "enum" in spec:
                options["choices"] = [value for value in spec["enum"] if value is not None]
            if name == 'groups':
                options['help'] = 'JSON array of {values: [...], replacement: "label"} rules; mutually exclusive with --values.'
            elif name == "values":
                options["help"] = 'JSON array (e.g. [1,2,null]); omitted or [] keeps all rows.'
            elif name == "conversion_rules":
                options["help"] = "JSON file containing [{from, to, formula}, ...]."
            elif "default" in spec:
                options["help"] = spec.get("description", f"Default: {spec['default']}.")
            else:
                options["help"] = "Required."
            command.add_argument("--" + name.replace("_", "-"), **options)
    return parser


def _json_load(text):
    def constant(value):
        raise ValueError(f"Non-finite JSON number: {value}")

    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f"Duplicate JSON key: {key}")
            result[key] = value
        return result
    try:
        return json.loads(text, parse_constant=constant, object_pairs_hook=pairs)
    except ValueError as exc:
        raise CliError("INVALID_JSON", str(exc)) from exc


def _input_spec(args):
    spec = {"path": str(Path(args.input).resolve()), "columns": {}}
    for name in ("header", "delimiter", "encoding", "sheet"):
        if hasattr(args, name):
            spec[name] = getattr(args, name)
    if hasattr(args, "sheet_index"):
        if args.sheet_index < 0:
            raise CliError("INVALID_ARGUMENTS", "--sheet-index must be nonnegative.")
        spec["sheet"] = args.sheet_index
    for role in COLUMN_ROLES:
        if hasattr(args, role + "_col"):
            spec["columns"][role] = getattr(args, role + "_col")
    if Path(spec["path"]).suffix.lower() not in INPUT_FORMATS:
        raise CliError("UNSUPPORTED_FORMAT", "Unsupported input suffix.", {"supported": INPUT_FORMATS})
    return spec


def _execute(args):
    if args.command == "schema":
        from .cli_runner import envelope
        result = envelope("schema")
        result["schema"] = schema_document()
        result["schema"]["exit_codes"] = {"0": "success", "1": "unexpected error", "2": "invalid arguments/configuration/data", "3": "file or dependency error", "4": "network failure or deadline exceeded", "5": "strict record failure", "130": "cancelled"}
        result["schema"]["response_format"] = "One JSON object on stdout; logs on stderr. --format text selects a human-readable response."
        return result
    protected = []
    rule_changes = None
    if getattr(args, "rules_file", None):
        path = Path(args.rules_file).resolve()
        rule_changes = normalize_rules(_json_load(path.read_text(encoding="utf-8-sig")))
        protected.append(path)
    if args.command == "rules":
        from .chem_processor import ChemistryProcessor
        from .cli_operations import apply_rules
        from .cli_runner import envelope
        effective = apply_rules(ChemistryProcessor(), rule_changes or {})
        return envelope("rules", rules=effective, summary={"counts": {key: len(value) for key, value in effective.items()}})
    if args.command == "column-values":
        from .cli_runner import envelope, json_safe, read_frame
        from .value_filter import column_values
        frame = read_frame(_input_spec(args))
        values = column_values(frame, args.column)
        return envelope("column-values", summary={"column": args.column, "rows": len(frame),
                                                  "distinct_values": len(values), "values": json_safe(values)})
    if args.command == "inspect":
        if not 0 <= args.limit <= 1000:
            raise CliError("INVALID_ARGUMENTS", "--limit must be between 0 and 1000.")
        from .cli_runner import inspect_file
        return inspect_file(_input_spec(args), args.limit)
    if args.command == "run":
        if args.config == "-":
            raw = _json_load(sys.stdin.read())
            base = Path.cwd()
        else:
            path = Path(args.config).resolve()
            raw = _json_load(path.read_text(encoding="utf-8-sig"))
            base = path.parent
            protected.append(path)
        if isinstance(raw, dict):
            raw = dict(raw)
            for flag in ("overwrite", "strict"):
                if getattr(args, flag, False):
                    raw[flag] = True
    else:
        params = {name: getattr(args, name) for name in PARAM_SPECS[args.command] if hasattr(args, name)}
        if args.command == 'transform-column':
            params.pop('unit_col', None)
            if hasattr(args, 'column_unit'):
                params['unit_col'] = args.column_unit
        if "values" in params:
            params["values"] = _json_load(params["values"])
        if 'groups' in params:
            params['groups'] = _json_load(params['groups'])
        if "conversion_rules" in params:
            rules_path = Path(params["conversion_rules"]).resolve()
            params["conversion_rules"] = _json_load(rules_path.read_text(encoding="utf-8-sig"))
            protected.append(rules_path)
        output = {"path": args.output}
        for option, key in (("report", "report_path"), ("excluded", "excluded_path"), ("columns", "columns")):
            if hasattr(args, option):
                output[key] = getattr(args, option)
        raw = {"schema_version": "1", "input": _input_spec(args), "steps": [{"op": args.command, "params": params}],
               "output": output, "overwrite": getattr(args, "overwrite", False), "strict": getattr(args, "strict", False)}
        base = Path.cwd()
    if rule_changes is not None and isinstance(raw, dict):
        raw["rules"] = rule_changes
    config = normalize_config(raw, base)
    from .cli_runner import execute
    return execute(config, args.command, getattr(args, "dry_run", False), protected)


def main(argv=None):
    """Return a process exit code; never request interactive input."""
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8")
    args = None
    exit_code = 0
    try:
        args = build_parser().parse_args(argv)
        # Importing core is deferred until after parsing and configuration validation.
        if args.command != "schema":
            from .logger import log_manager
            log_manager.configure(enable_file=False, console_level=logging.INFO if args.verbose else logging.WARNING)
        with redirect_stdout(sys.stderr):
            result = _execute(args)
    except (Exception, KeyboardInterrupt) as exc:
        if args is not None and args.debug:
            import traceback
            traceback.print_exc(file=sys.stderr)
        if not isinstance(exc, CliError):
            if isinstance(exc, KeyboardInterrupt):
                exc = CliError("CANCELLED", "Operation cancelled; no completion report was published.", exit_code=130)
            elif isinstance(exc, ModuleNotFoundError):
                exc = CliError("MISSING_DEPENDENCY", str(exc), {"module": exc.name}, exit_code=3)
            elif isinstance(exc, OSError):
                exc = CliError("IO_ERROR", str(exc), exit_code=3)
            elif isinstance(exc, (ValueError, KeyError, UnicodeError, LookupError)):
                exc = CliError("INVALID_DATA", str(exc))
            else:
                exc = CliError("INTERNAL_ERROR", str(exc), exit_code=1)
        exit_code = exc.exit_code
        from .cli_runner import envelope
        result = envelope(args.command if args else None, status="error", error={
            "code": exc.code, "message": exc.message, "details": exc.details, "retryable": exc.retryable,
        })
    if args is not None and args.format == "text":
        if result["status"] == "error":
            print(f"{result['error']['code']}: {result['error']['message']}")
            if result["error"]["details"]:
                print(json.dumps(result["error"]["details"], ensure_ascii=False, indent=2))
        else:
            print(f"{result['command']}: success")
            print(json.dumps(result.get("schema", result.get("rules", result["summary"])), ensure_ascii=False, allow_nan=False, indent=2))
            for key, path in result["artifacts"].items():
                print(f"{key}: {path}")
            for warning in result["warnings"]:
                print(f"Warning: {json.dumps(warning, ensure_ascii=False)}")
    else:
        print(json.dumps(result, ensure_ascii=False, allow_nan=False))
    return exit_code


if __name__ == "__main__":
    from multiprocessing import freeze_support
    freeze_support()
    raise SystemExit(main())
