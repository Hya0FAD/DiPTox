# DiPTox - Data Integration and Processing for Computational Toxicology

[![PyPI](https://img.shields.io/pypi/v/diptox)](https://pypi.org/project/diptox/) [![Conda](https://img.shields.io/conda/vn/conda-forge/diptox.svg)](https://anaconda.org/conda-forge/diptox) [![Conda Platforms](https://img.shields.io/conda/pn/conda-forge/diptox.svg)](https://anaconda.org/conda-forge/diptox) ![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg) ![Python Version](https://img.shields.io/badge/python-3.10+-brightgreen.svg) [![Chinese](https://img.shields.io/badge/-%E4%B8%AD%E6%96%87%E7%89%88-blue.svg)](./README_ZH.md) [![PyPI Downloads](https://static.pepy.tech/badge/diptox)](https://pepy.tech/project/diptox) [![Conda Downloads](https://img.shields.io/conda/dn/conda-forge/diptox.svg)](https://anaconda.org/conda-forge/diptox)
<p align="center">
  <img src="assets/TOC.png" alt="DiPTox Workflow Diagram" width="500">
</p>
**DiPTox** is a Python toolkit designed for the robust preprocessing, standardization, and multi-source data integration of molecular datasets, with a focus on computational toxicology workflows.

## v1.1.0 Updates

- **CLI and JSON pipelines**: Added a command-line entry point for agents and scripts, with discovery, data inspection, configuration preflight, preprocessing, unit conversion, deduplication, substructure search, atom-count filtering, and InChI calculation.
- **Network enrichment and chemical rules**: Query multiple sources with rate limits, retries, a total deadline, and field provenance; configure chemical rules for each invocation through JSON.
- **NiceGUI interface**: Keep settings across pages and run long operations in the background. Jobs that modify data commit their results only after completion.
- **Results and audit**: The CLI provides structured JSON responses, defined exit codes, source-row tracking, exclusion records, and completion reports for reproducible processing.

See the [CLI Guide](docs/CLI.md) for usage and the [Changelog](CHANGELOG.md#english) for current and previous version notes.

## DiPTox Community Check-in (Optional)
To help us understand our user base and improve the software, DiPTox includes a one-time, optional survey on first use. 
-   **Completely Optional**: You can skip it with a single click.
-   **Privacy-Focused**: The information helps us with academic impact assessment. It will not be shared.

## Core Features

#### 1. Graphical User Interface (GUI)
Powered by NiceGUI, the local web interface allows users to perform all workflows visually without writing code.
-   **Persistent Pages**: Configuration values remain intact while navigating between workflow steps.
-   **Responsive Jobs**: Long-running processing and network operations do not block page navigation.
-   **Real-time Preview**: Instantly view data changes after applying rules.
-   **Rule Management**: Add/Remove valid atoms, salts, solvents, and **unit conversion formulas** interactively.
-   **Smart Column Mapping**: Intelligent detection of headers and binary file structures.

#### 2. Chemical Preprocessing & Standardization
A configurable pipeline to clean and normalize chemical structures.
-   **Strict Inorganic Filtering**: Updated SMARTS matching to accurately identify complex inorganic species (e.g., ionic cyanides) without misclassifying organic nitriles.
-   **Pipeline Steps**:
    -   Remove salts & solvents
    -   Handle mixtures (keep largest fragment)
    -   Remove inorganic molecules
    -   Neutralize charges & Validate atomic composition
    -   Remove explicit hydrogens, stereochemistry, and isotopes
    -   **Reject Radical Species**: Automatically discard molecules containing free radical atoms.
    -   Standardize to canonical SMILES
    -   Filter by atom count

#### 3. Unit Standardization
Normalize heterogeneous target data into a single unit effortlessly.
-   **Automatic Conversion**: Built-in rules for **Concentration**, **Time**, **Pressure**, and **Temperature**.
-   **Custom Formulas**: Define mathematical rules (e.g., `x * 1000` or `10**(-x)`) interactively via GUI or script.
-   **Unified Output**: Standardize diverse units (e.g., `ug/mL`, `g/L`, `M`) to a single target (e.g., `mg/L`).

#### 4. Data Deduplication
Flexible strategies for handling duplicate entries with advanced controls.
-   **Data Types**: Supports `continuous` (e.g., IC50) and `discrete` (e.g., Active/Inactive) targets.
-   **Continuous values**: `method="auto"`, `"IQR"`, or `"3sigma"` selects outlier filtering; `aggregation="mean"`, `"max"`, or `"min"` selects aggregation of the remaining values. Defaults: auto + mean. Groups of at most 3 skip filtering.
-   **Discrete values**: `method="vote"` or `"priority"`.
-   **Log Transformation**: Optional `-log10` transformation (e.g., IC50 $\to$ pIC50) applied *before* deduplication logic to handle bioactivity data correctly.
-   **Flexible NaN Handling**: Option to retain rows with missing conditions (treating *NaN* as a valid group) instead of dropping them.

#### 5. Comprehensive History Tracking (Audit Log)
-   Records every operation (Loading, Preprocessing, Filtering, etc.) in an **Audit Log**.
-   Tracks **timestamps**, **operation details**, and row count changes (**Delta**) step-by-step.
-   Available via API (`get_history()`) and visualized in the GUI.

#### 6. Identifier & Property Integration
-   Fetch and interconvert identifiers (**CAS, SMILES, IUPAC, MW**) from multiple sources (**PubChem, ChemSpider, CompTox, Cactus, CAS Common Chemistry, ChEMBL**).
-   High-performance **concurrent requests** with automatic rate limiting and retries.

#### 7. Utility Tools
-   Perform **substructure searches** using SMILES or SMARTS patterns.
-   **Customize chemical processing rules** for neutralization reactions, salt/solvent lists, and valid atoms.
-   **Display a summary** of all currently active processing rules.

## Installation
You can install DiPTox using `pip` or via `conda`/`mamba`.
### Option 1: Install from PyPI
Install the official stable version from PyPI:
```bash
pip install diptox
```
### Option 2: Install from Conda-forge
Installing `diptox` from the `conda-forge` channel can be achieved by adding `conda-forge` to your channels with:

```bash
conda config --add channels conda-forge
conda config --set channel_priority strict
```

Once the `conda-forge` channel has been enabled, `diptox` can be installed with `conda`:

```bash
conda install diptox
```

or with `mamba`:

```bash
mamba install diptox
```

## GUI
After installation, you can launch the graphical interface directly from your terminal:

```bash
diptox-gui
```

This command starts DiPTox locally and opens the interface in your default browser.

## Command-line interface (CLI)

Use the CLI for agents, scripts, and reproducible JSON pipelines. Start with `diptox --help` or `python -m diptox --help`.

See the [CLI Guide](docs/CLI.md) ([Chinese](docs/CLI_ZH.md)) for commands, configuration, chemical rules, network enrichment, and audit outputs.

## Quick Start
```python
from diptox import DiptoxPipeline

def main():
    # Initialize processor
    DP = DiptoxPipeline()

    # Load data
    DP.load_data(input_data='file_path/list/dataframe', smiles_col, target_col, cas_col, unit_col)

    # Customize Processing Rules (Optional)
    print("--- Default Rules ---")
    DP.display_processing_rules()

    DP.manage_atom_rules(atoms=['Si'], add=True)         # Add 'Si' to the list of valid atoms
    DP.manage_default_salt(salts=['[Na+]'], add=False)   # Example: remove sodium from the salt list
    DP.manage_default_solvent(solvents='Cl', add=False)  # Example: remove chlorine from the solvent list
    DP.add_neutralization_rule('[$([N-]C=O)]', 'N')      # Add a custom neutralization rule

    print("\n--- Customized Rules ---")
    DP.display_processing_rules()

    # Configure preprocessing
    DP.preprocess(
        remove_salts=True,              # Remove salt fragments. Default: True.
        remove_solvents=True,           # Remove solvent fragments. Default: True.
        mixture_mode="reject",          # keep / reject / largest. Default: reject.
        hac_threshold=3,                # Largest mode requires a unique largest fragment with HAC > 3.
        remove_inorganic=False,         # Remove common inorganic molecules. Default: True.
        neutralize=True,                # Neutralize charges on the molecule. Default: True.
        reject_non_neutral=False,       # Only retain the molecules whose formal charge is zero. Default: False.
        element_policy="allow_all",     # allow_all / reject_metals / allowed_atoms. Default: allow_all.
                                        # allowed_atoms rejects the whole molecule if any element is unlisted, regardless of atom degree.
        remove_stereo=False,            # Remove stereochemistry information. Default: False.
        remove_isotopes=True,           # Remove isotopic information. Default: True.
        remove_hs=True,                 # Remove explicit hydrogen atoms. Default: True.
        add_hs=False,                   # Add explicit H after all chemical processing. Default: False.
        reject_radical_species=True,    # Molecules containing free radical atoms are directly rejected. Default: True.
        n_jobs=4                        # Accelerate using 4 CPU cores. Default: 1 
    )

    # Configure deduplication and unit standardization
    conversion_rules = {('g/L', 'mg/L'): 'x * 1000', 
                        ('M', 'mg/L'): 'x * mw * 1000',}
    DP.standardize_units(standard_unit="mg/L", conversion_rules=conversion_rules,
                         molecular_weight_source="original")
    DP.config_deduplicator(condition_cols=condition_cols, data_type=data_type, method=method, aggregation="mean")
    DP.dataset_deduplicate()

    # Configure web queries
    DP.config_web_request(sources=['pubchem/chemspider/comptox/cactus/cas'], max_workers, ...)
    DP.web_request(send='cas', request=['smiles', 'iupac'])

    # Substructure search
    DP.substructure_search(query_pattern, is_smarts=True)

    # Save results
    DP.save_results(output_path='file_path')

    # View Processing History (Audit Log)
    print(DP.get_history())
    # Output Example:
    #               Step Timestamp  Rows Before  Rows After   Delta                               Details
    # 0     Data Loading  10:00:01            0        1000   +1000                   Source: dataset.csv
    # 1    Preprocessing  10:00:05         1000         950     -50  Valid: 950, Invalid: 50. Order: ...
    # 2    Deduplication  10:00:08          950         800    -150       Method: auto (Log10 Transformed)

# CRITICAL: This protection block is REQUIRED for Windows multiprocessing!
# It prevents infinite recursive loops and memory explosion when n_jobs > 1.
if __name__ == '__main__':
    main()
```

Choose `molecular_weight_source="original"` or `"standardized"`, or provide `molecular_weight_col`, whenever a formula uses `mw`. `standardized` requires preprocessing first; the source remains unset by default. Unit scaling that does not use MW needs no selection. `mixture_mode="largest"` rejects tied largest fragments as `Ambiguous parent`. When `remove_hs` and `add_hs` are both enabled, hydrogen removal happens first and addition happens after all other chemical processing. The Python API retains older mixture/element switches for compatibility; new code should use `mixture_mode` and `element_policy`.

`Canonical SMILES` omits atom-map numbers; the input SMILES and `Original Canonical SMILES` retain source annotations. Salt removal, solvent removal and mixture handling each collapse repeated components only when all remaining components have the same identity: `A.A` becomes `A`, while `A.A.B` is not reduced to `A.B`. Mixture modes apply after this normalization, including `keep`. There is no separate switch. Repeated copies of a single known solvent also collapse to one copy; other selected chemistry rules still apply.

Ethylene glycol (`OCCO`) and 2-methoxyethanol (`COCCO`) are solvent rules, not salt rules. Salt and solvent recognition share whole-fragment matching with support for the `S=O` / `[S+][O-]` representations of sulfoxides such as DMSO. `reject_metals` uses an explicit metal set; noble gases and B, Si, Ge, As, Sb and Te are outside that set. `allowed_atoms` separately checks the user-configured atom list.

## Advanced Configuration

### Repeated processing and source data

When an imported field conflicts with a generated column, DiPTox preserves it as `Original name (Input)` and updates its mapping automatically. No spreadsheet renaming is needed. SDF export synchronizes explicitly mapped custom structure fields and archives their original declarations.

Neutralization chooses sites in canonical order so equivalent SMILES yield the same result. When component processing is enabled, identical components are collapsed again after neutralization or hydrogen addition; `A.A.B` is still not partially deduplicated. Non-string SMILES are rejected per record without interrupting the batch.

If reprocessing changes the standardized molecular weight, targets depending on the old MW are cleared and marked `Stale standardized molecular weight`. Reconvert from the original concentrations before deduplication; aggregated data requires undo or restoration of the original records. Pure scaling and MW sources independent of the standardized structure remain valid.

`log_transform` preserves input values and writes a separate target with units such as `-log10(M)`. Repeating the same deduplication does not apply the logarithm again. Linear unit conversion rejects logarithmic targets; undo the transformation or select the original linear data first. Newly requested deduplication units are applied, and continuous targets containing `inf/-inf` are excluded.

### Web Service Integration
DiPTox supports the following chemical databases:
-   `PubChem`: https://pubchem.ncbi.nlm.nih.gov/
-   `ChemSpider`: https://www.chemspider.com/
-   `CompTox`: https://comptox.epa.gov/dashboard/
-   `Cactus`: https://cactus.nci.nih.gov/
-   `CAS`: https://commonchemistry.cas.org/
-   `ChEMBL`: https://www.ebi.ac.uk/chembl/

**Note:** `ChemSpider`, `CompTox` and `CAS` require API keys. Provide them during configuration:
```python
DP.config_web_request(
    sources=['chemspider/comptox/CAS'],
    chemspider_api_key='your_personal_key',
    comptox_api_key='your_personal_key',
    cas_api_key='your_personal_key'
)
```
## Requirements
- `Python>=3.8`
- Core Dependencies:
  - `requests`
  - `rdkit>=2023.3`
  - `tqdm`
  - `openpyxl`
  - `scipy`
  - Python 3.8/3.9: `nicegui==2.24.2`; Python 3.10+: `nicegui>=3.16,<4` (selected automatically during installation)
- Optional Dependencies (install as needed, if not installed, then send the request using `requests`.):
  - `pubchempy>=1.0.5`: For PubChem integration
  - `chemspipy>=2.0.0`: For ChemSpider (requires API key)
  - `ctx-python>=0.0.1a10`: For CompTox Dashboard (requires API key)

Both Python branches share the current GUI and processing features, including background jobs, chemical rules, unit conversion, deduplication, filtering, and export. Compatibility adapters handle legacy uploads and Python 3.8 background thread execution. Run `python -m pip install .` in the desired Python environment to install this source checkout. Python 3.8/3.9 uses the legacy framework; the GUI binds to localhost (`127.0.0.1`) by default.

## License
Apache License 2.0 - See [LICENSE](LICENSE) for details

## Support
Report issues on [GitHub Issues](https://github.com/Hya0FAD/DiPTox/issues)
