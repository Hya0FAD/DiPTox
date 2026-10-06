# DiPTox CLI Guide

[Back to README](../README.md) | [Chinese guide](CLI_ZH.md)

Unless stated otherwise, commands assume the repository root is the working directory; example file paths are relative to that root.

The CLI is designed for agents and unattended scripts. Run `diptox` after installing this version, or use `python -m diptox` from the repository. For a source checkout, `python -m pip install -e .` also makes it available from other directories. `diptox-gui` still opens the GUI. Commands cover discovery, inspection, preprocessing, unit conversion, deduplication, substructure search, atom-count filtering, InChI calculation, chemical rules, and bounded network enrichment. Combine processing commands as ordered steps in a JSON pipeline.

## Inspect, configure, execute

From the repository root:

```bash
python -m diptox schema
python -m diptox inspect --input examples/cli/data.csv --limit 3
python -m diptox run --config examples/cli/pipeline.json --dry-run
python -m diptox run --config examples/cli/pipeline.json
```

The [example configuration](../examples/cli/pipeline.json) processes [six offline records](../examples/cli/data.csv), including two duplicate pairs, mixed units, an invalid structure, and an invalid numeric value. It writes `clean.csv`, `excluded.csv`, and `report.json` under `examples/cli/output/`. The final dataset has two rows, with standardized values of 1 and 2 mg/L, and the audit contains four exclusion events. Add `--overwrite` to rerun after reviewing the outputs.

`schema` returns the configuration schema, step parameters and defaults; `inspect` reports columns, dtypes and a bounded preview. Inspection loads the selected table; `--limit` limits the displayed preview, not input reading. For workbooks with several sheets, inspection first lists their names; select one with `--sheet "Sheet1"` or the zero-based `--sheet-index 0`.

`--dry-run` validates configuration, input metadata, mapped and generated columns, chemical rules, query syntax, unit-rule availability and output conflicts. For enrichment it also checks provider capabilities and required environment credentials, without making network requests. It reads the input but does not run the processing steps, validate every molecule, guarantee successful conversions or remote results, or write artifacts.

## Individual steps and configuration

```bash
diptox preprocess --input data.csv --smiles-col SMILES --output processed.csv --excluded rejected.csv --report preprocess.json
diptox units --input data.csv --target-col Value --unit-col Unit --standard-unit mg/L --output units.csv
diptox deduplicate --input processed.csv --smiles-col "Canonical SMILES" --data-type smiles --output unique.csv
```

Use `diptox COMMAND --help` for all flags. Column roles are `smiles`, `target`, `unit`, `cas`, `name`, `inchikey`, and `id`, with corresponding `--smiles-col`, `--target-col`, etc. Map the columns required by each operation explicitly. Headerless input uses `--no-header` and column names `0`, `1`, etc.; for example, `--smiles-col 0 --target-col 1 --unit-col 2`. Text input defaults to UTF-8; CSV uses a comma and TXT a tab unless `--delimiter` is supplied.

JSON configurations require `"schema_version": "1"`, `input`, a nonempty ordered `steps` array, and `output`. A step has `op` and `params`; input mappings use `input.columns`, and artifact paths use `output.path`, optional `excluded_path`, and optional `report_path`. Unknown keys and incompatible settings fail validation. Defaults include:

| Setting | Default / supported values |
| --- | --- |
| Input | `header: true`, `encoding: "utf-8"` |
| Execution | `overwrite: false`, `strict: false`; no interactive prompts |
| Preprocessing | `mixture_mode: "reject"`, `element_policy: "allow_all"`, `n_jobs: 1`, `chunksize: 100`, `hac_threshold: 3` |
| Chemical flags | Salt/solvent/inorganic removal, neutralization, isotope/H removal, sanitization and radical rejection enabled; stereo removal, H addition (`add_hs`), non-neutral rejection and strict atom checking disabled |
| Units | `standard_unit` required; molecular-weight basis initially unset; choose `molecular_weight_source` or `molecular_weight_col` when a conversion uses MW |
| Deduplication | `data_type: "continuous"`, `method: "auto"`, `p_threshold: 0.05`, `log_transform: "None"`, `condition_cols: []`, `dropna_conditions: false` |

Boolean step options support both forms, such as `--remove-salts` / `--no-remove-salts`. Mixtures support `keep`, `reject`, `largest`; element policies support `allow_all`, `reject_metals`, `allowed_atoms`. Continuous deduplication uses `--method auto|IQR|3sigma` for outlier filtering and `--aggregation mean|max|min` for aggregation after filtering (JSON keys: `method`, `aggregation`). Defaults are auto + mean; groups of at most 3 skip filtering. Old `method: "max"` / `"min"` configurations must move the extremum choice to `aggregation`; discrete data defaults to `vote` and also supports `priority` with a nonempty ordered list. Structure-only deduplication uses `data_type: "smiles"` and `method: "auto"`. Log transforms accept the strings `"None"`, `"-log10"`, `"log10"` and apply only to continuous data; pass a leading-minus CLI value as `--log-transform=-log10`.

Custom unit rules use `"conversion_rules": [{"from": "custom", "to": "mg/L", "formula": "x * 1000"}]`; the standalone command accepts the same array from `--conversion-rules rules.json`. Rules use `x`, optional molecular weight `mw`, arithmetic, and `log`, `log10`, `exp`. Every conversion using `mw`, including mass/molar conversions, requires an explicit `molecular_weight_source` of `original` or `standardized`, or a `molecular_weight_col`. `standardized` requires an earlier preprocessing step. The `auto` MW option is no longer accepted. Pure unit scaling can leave the basis unset.

`add_hs: true` (CLI `--add-hs`) adds explicit hydrogens after all chemical processing. It can be combined with `remove_hs: true`: removal happens first, addition happens last. In `mixture_mode: "largest"`, tied largest fragments are rejected as `Ambiguous parent`; no fragment is chosen arbitrarily. Prefer `mixture_mode` and `element_policy` in Python as well; the Python API retains the older switches only for compatibility.

Preprocessing always removes atom-map numbers from `Canonical SMILES`, retaining the input annotations in the original columns. Map-only differences do not count as different parent structures or a chemical structure change. Salt removal, solvent removal and mixture handling all collapse an all-identical set of remaining components (`A.A` → `A`); none performs partial duplicate removal (`A.A.B` → `A.B`). The identity comparison ignores atom maps but retains the stereochemistry, isotopes, charge and hydrogen representation present at that step. Repeats consisting entirely of one known solvent also collapse. Mixture normalization precedes all three modes (`keep`, `reject`, `largest`); a repeated single identity is not subject to the largest-fragment size threshold. Mixed identities still follow the chosen mixture policy. No new flag is required.

Ethylene glycol (`OCCO`) and 2-methoxyethanol (`COCCO`) are controlled by `remove_solvents`, not by salt removal. Standalone molecules are retained. Salt and solvent recognition share whole-fragment matching that accepts both `S=O` and `[S+][O-]` forms of ordinary sulfoxides with two carbon substituents. These transformations affect matching copies only; atom counts, bond counts and query constraints still apply to the entire fragment. Larger molecules containing the same functional group are not removed as DMSO, and no global tautomer or charge normalization is performed.

`reject_metals` uses an explicit set of alkali, alkaline-earth, d-block, lanthanide, actinide and listed p-block metals, independent of the allowed-atom list. Noble gases, halogens and B/Si/Ge/As/Sb/Te are outside the metal set. The boundary element Po is classified as a metal by this application's convention; see `diptox/element_policy.py` for the classification and sources. Rejection reasons, salt-fragment protection and `Original/Final Metal Elements` share this definition.

Paths inside a configuration file are relative to that file's directory. `run --config -` reads JSON from stdin, with paths relative to the working directory. For example, from `examples/cli`, use `cat pipeline.json | python -m diptox run --config -` on POSIX shells, or `Get-Content -Raw pipeline.json | python -m diptox run --config -` in PowerShell. Standalone command paths are also relative to the working directory. `run --overwrite` and `run --strict` enable the corresponding configuration settings.

## Search, filtering and identifiers

```bash
diptox search --input data.csv --smiles-col SMILES --query-pattern "[OX2H]" --query-type smarts --mode matches --output hydroxyl.csv
diptox search --input data.csv --smiles-col SMILES --query-pattern "CCO" --query-type smiles --mode annotate --output annotated.csv
diptox filter-atoms --input data.csv --smiles-col SMILES --min-heavy-atoms 3 --max-heavy-atoms 30 --max-total-atoms 100 --output filtered.csv
diptox inchi --input data.csv --smiles-col SMILES --output identifiers.csv
```

`search` requires `query_pattern`; `query_type` defaults to `smarts` and also accepts `smiles`. All modes add `Substructure_<query_pattern>`: `annotate` (default) retains every row, `matches` keeps matching structures, and `nonmatches` keeps nonmatching structures. Invalid or rejected structures have an unknown (empty) result and are excluded in either selection mode; they are not counted as nonmatches.

`filter-atoms` accepts `min_heavy_atoms`, `max_heavy_atoms`, `min_total_atoms` and `max_total_atoms`. Bounds are inclusive nonnegative integers; supply at least one bound, with each minimum no greater than its maximum. Total atom counts include hydrogen atoms. `inchi` adds an `InChI` column and records invalid structures or generation failures in the audit. After preprocessing, these operations use the standardized structures.

The [offline P2 example](../examples/cli/pipeline-p2.json) reuses `data.csv`: preprocessing → oxygen substructure selection → three-heavy-atom filtering → InChI → structure-only deduplication. It produces one ethanol record under `examples/cli/output-p2/` and audits invalid structures and ordinary filtering separately.

```bash
python -m diptox run --config examples/cli/pipeline-p2.json --dry-run
python -m diptox run --config examples/cli/pipeline-p2.json
```

## Chemical rules for one invocation

```bash
diptox rules
diptox rules --rules-file examples/cli/chemical-rules.json
diptox preprocess --input data.csv --smiles-col SMILES --element-policy allowed_atoms --rules-file examples/cli/chemical-rules.json --output processed.csv
diptox run --config examples/cli/pipeline-p2.json --rules-file examples/cli/chemical-rules.json
```

`rules` returns the effective chemical rules as JSON. Every processing command and `run` accepts `--rules-file`; the file contains rule changes directly, while pipeline JSON stores the same object in a top-level `rules` field. A command-line rules file **replaces** that field rather than merging with it. Changes affect only the current invocation, without changing global defaults or creating persistent settings.

```json
{
  "atoms": {"add": ["Xe"], "remove": []},
  "salts": {"add": ["[Na+]"], "remove": []},
  "solvents": {"add": ["CCOCC"], "remove": []},
  "neutralization": {
    "add": [{"reactant": "[O-;X1]", "product": "O"}],
    "remove": []
  }
}
```

Sections and change lists are optional. Atoms use element symbols; salts use SMARTS; solvents use SMILES; neutralization additions pair reactant SMARTS with replacement SMILES, and removals identify a reactant SMARTS rule. Removals are applied before additions; removing an absent rule fails explicitly. Preflight validates chemical syntax. Atom-list changes take effect when preprocessing uses `element_policy: "allowed_atoms"`. The report includes effective rule counts and a SHA-256 fingerprint for reproducibility.

## Network enrichment

```bash
diptox enrich --input examples/cli/data-enrich.csv --name-col Name --sources pubchem --send name --request smiles cas iupac mw --timeout 15 --deadline 120 --output enriched.csv --report enrich-report.json
python -m diptox run --config examples/cli/pipeline-enrich.json --dry-run
```

The [network example](../examples/cli/pipeline-enrich.json) starts with names, requests structures and properties, then preprocesses the returned SMILES and calculates InChI. Run it without `--dry-run` to contact the selected service. It writes to `examples/cli/output-enrich/`. Input mappings describe the original file: when an enrichment step requests `smiles`, `cas` or `name`, subsequent steps can use the generated role. A SMILES request after preprocessing preserves the existing standardized structure role.

| Parameter | Default / supported values |
| --- | --- |
| `sources` | `["pubchem"]`; ordered fallback across `pubchem`, `chemspider`, `comptox`, `cactus`, `chembl`, `cas` |
| `send` | `["smiles"]`; ordered identifiers from `smiles`, `cas`, `name`; map every selected role |
| `request` | Required nonempty list from `smiles`, `cas`, `iupac`, `mw`, `name` |
| `max_workers` | `4`, from `1` through `32` |
| `timeout` | HTTP transport timeout, default `15` seconds per request |
| `deadline` | `120` seconds for the complete enrichment worker, including startup, requests, retries and waits |
| `retries` | `2` additional attempts for transient HTTP failures; `0` sends only the initial attempt |
| `retry_delay` | `1` second between retry attempts |
| `interval` | `0.3` seconds minimum between HTTP request starts across the worker |

CLI enrichment always uses direct HTTP APIs, including when optional provider SDKs are installed. Unsupported provider/property or identifier combinations fail preflight. ChemSpider, CompTox and CAS credentials are read from `DIPTOX_CHEMSPIDER_API_KEY`, `DIPTOX_COMPTOX_API_KEY` and `DIPTOX_CAS_API_KEY`, respectively. Set the relevant environment variable before selecting that provider; keys are not accepted in configuration files or written to reports. For example, in PowerShell:

```powershell
$env:DIPTOX_CHEMSPIDER_API_KEY = "your-key"
```

Returned properties use `<property>_from_web` columns. `Query_Status` records `complete`, `partial`, `not_found`, `invalid_identifier` or `failed`; `Query_Missing_Fields` and `Query_Errors` contain JSON arrays. `Data_Source` contains a JSON object keyed by property, recording its provider and identifier type; `Query_Method` records successful identifier types. Step summaries count `requested_fields`, `resolved_fields` and `missing_fields` across all row/property cells, plus `complete_rows`, `partial_rows`, `failed_rows` and classified `error_counts`.

Missing or partial results produce row-level warnings and audit events. If no fields are resolved and network/provider failures occurred, the command returns exit code `4` without publishing artifacts. Exceeding `deadline` terminates the worker and its request threads, returns `DEADLINE_EXCEEDED` with exit code `4`, and publishes no artifacts. Cancellation returns `130` and stops the worker. `--strict` turns incomplete row results into exit code `5` before publication; ordinary `not_found` results do not by themselves imply a network outage.

## Results, audit and failure handling

`run` keeps column/state changes between steps. Units create `Value (Standardized)` and `Unit (Standardized)` in this example; deduplication creates `Value (Standardized)_new`. A CSV does not restore the pipeline's internal state when reloaded. For separate commands, explicitly map the actual generated names, including `--smiles-col "Canonical SMILES"` and the current target/unit columns. Use the report's `summary.column_roles` to discover the final names.

Preprocessing retains rejected rows in the main table, marked by `Is Valid`, `Standardization Status`, and `Processing Log`. Unit failures also remain with status columns; deduplication removes unusable rows and consolidates duplicate groups. The reserved `Diptox Source Rows` column tracks zero-based input record positions (excluding a header), joined with semicolons after consolidation. Exclusion exports contain per-step events with `CLI Step`, `CLI Operation`, `CLI Reason` and `CLI Severity`: one source record may appear at several steps, so `exclusion_events` is not a count of unique rejected records. Normal search/atom-count exclusions have severity `filtered`; they contribute to `filtered_rows` and do not fail strict mode. Invalid structures have severity `error` and still fail strict mode. Consolidated duplicates are reported separately from failures.

Record failures normally produce warnings and audited output with exit code 0. `--strict` makes these failures exit with code 5 before publishing any artifacts; existing files remain untouched. Inputs, configuration files and outputs must have distinct paths. Outputs are not overwritten without `--overwrite`. All artifacts are staged first, then published atomically **per file**; this is not a multi-file transaction. The completion report is published last. When overwriting, an old completion report is set aside before publication and is not retained if only part of the new run is committed. A publication failure can leave earlier outputs committed; inspect `error.details.published_artifacts` and the process result before treating a run as complete.

Inputs: CSV, XLSX, XLS, TXT, SMI, SDF, MOL. Main outputs: CSV, XLSX, TXT, SMI, SDF; XLS export is not supported. XLS input needs a compatible pandas Excel engine installed. Exclusions use CSV/XLSX/TXT; reports use JSON. Python, GUI and CLI use the same SDF loader: failed records remain as invalid rows with a one-based `Input Record` and `Import Error`. SDF output preserves invalid rows with an empty molecule block and `SDF Export Status: Invalid structure placeholder`; their data and audit properties remain available. SMI export requires valid structures. Optional `--columns` / `output.columns` selects result columns, including `Original Structure Type` after preprocessing.

The SDF molblock is the authoritative structure. A conflicting or invalid nonempty SMILES declaration is marked `Structure Reliability=Unreliable` and removed from the working dataset at load time, with its original declaration and reason retained in the exclusion audit. CLI reports these as step 0 (`load`); strict mode blocks publication. Comparison ignores atom-map numbers and ordinary explicit hydrogens but preserves stereochemistry, isotopes, charge and component multiplicity. Export rebuilds the molecule from the selected structure column, synchronizes generic SMILES properties and retains their original values in `SDF Source SMILES`.

Pipeline rows have stable private identities, so duplicate input index labels cannot mix results or exclusions. These identities are not output columns or GUI fields. Displayed source labels remain the original labels; repeated deduplication preserves the original membership. Direct use of the lower-level `DataDeduplicator` requires a unique DataFrame index; `DiptoxPipeline.load_data` assigns identities automatically.

Commands emit one UTF-8 JSON object on stdout with `schema_version`, `status`, `summary`, `artifacts`, `warnings`, and `error`. Failures include `error.code`, `message`, `details`, and `retryable`. Logs use stderr; `--verbose` enables informational logs and `--debug` adds failure tracebacks. `--format text` changes the response presentation, not the exported data format. `--help` and `--version` emit ordinary text.

| Exit code | Meaning |
| --- | --- |
| 0 | Success, possibly with record-level warnings |
| 1 | Unexpected internal error |
| 2 | Invalid arguments, configuration or data |
| 3 | File, publication or missing dependency error |
| 4 | Network/provider failure or enrichment deadline exceeded |
| 5 | Strict-mode record failure; no artifacts published |
| 130 | Cancelled |

**Logging compatibility:** importing DiPTox no longer creates a log directory or changes the application's root logger. File logging is opt-in, and the CLI disables it. Python applications that relied on automatic files should explicitly configure a `LogManager`, for example:

```python
from diptox import LogManager

LogManager().configure(enable_file=True)
```
