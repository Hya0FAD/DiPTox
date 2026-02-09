# DiPTox - Data Integration and Processing for Computational Toxicology

[![PyPI](https://img.shields.io/pypi/v/diptox)](https://pypi.org/project/diptox/) [![Conda Version](https://img.shields.io/conda/vn/conda-forge/diptox.svg)](https://anaconda.org/conda-forge/diptox) [![Conda Platforms](https://img.shields.io/conda/pn/conda-forge/diptox.svg)](https://anaconda.org/conda-forge/diptox) ![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg) ![Python Version](https://img.shields.io/badge/python-3.8+-brightgreen.svg) [![Chinese](https://img.shields.io/badge/-%E4%B8%AD%E6%96%87%E7%89%88-blue.svg)](./README_ZH.md)
<p align="center">
  <img src="assets/TOC.png" alt="DiPTox Workflow Diagram" width="500">
</p>
**DiPTox** is a Python toolkit designed for the robust preprocessing, standardization, and multi-source data integration of molecular datasets, with a focus on computational toxicology workflows.

## Official Release v1.0
We are excited to announce the first official stable release of DiPTox on PyPI! This milestone brings production-ready stability and performance enhancements:

* **Multi-Process Acceleration**:
    * Accelerate chemical preprocessing tasks by **10x or more** using the `n_jobs` parameter.
    * Intelligent task distribution across CPU cores for heavy datasets.
* **Cross-Platform Robustness**:
    * Implemented a specialized **"Guard Mechanism"** for Windows multiprocessing to prevent memory explosion and recursive process spawning issues.
    * Verified stability across Windows, Linux, and macOS environments.
* **Enhanced Data Loading**:
    * Switched to binary stream parsing for `.sdf` and `.mol` files to resolve encoding crashes (e.g., `utf-8` vs `latin-1`).
    * Auto-parsing of molecular structures to generate SMILES even when properties are missing.

## DiPTox Community Check-in (Optional)
To help us understand our user base and improve the software, DiPTox includes a one-time, optional survey on first use. 
-   **Completely Optional**: You can skip it with a single click.
-   **Privacy-Focused**: The information helps us with academic impact assessment. It will not be shared.

## Core Features

#### 1. Graphical User Interface (GUI)
Powered by Streamlit, the GUI allows users to perform all workflows visually without writing code.
-   **Visual Operation**: Complete workflow control via a web browser.
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
-   **Methods**: `auto`, `IQR`, `3sigma`, `vote`, or custom priority rules.
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

This command will automatically open the DiPTox interface in your default web browser.

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
        remove_mixtures=True,           # Handle mixtures based on fragment size. Default: False.
        hac_threshold=3,                # Heavy atom count threshold for fragment removal. Default: 3.
        keep_largest_fragment=True,     # Keep the largest fragment in a mixture. Default: True.
        remove_inorganic=False,         # Remove common inorganic molecules. Default: True.
        neutralize=True,                # Neutralize charges on the molecule. Default: True.
        reject_non_neutral=False,       # Only retain the molecules whose formal charge is zero. Default: False.
        check_valid_atoms=True,         # Check if all atoms are in the valid list. Default: False.
        strict_atom_check=False,        # If True, discard molecules with invalid atoms. If False, try to remove them. Default: False.
        remove_stereo=False,            # Remove stereochemistry information. Default: False.
        remove_isotopes=True,           # Remove isotopic information. Default: True.
        remove_hs=True,                 # Remove explicit hydrogen atoms. Default: True.
        reject_radical_species=True,    # Molecules containing free radical atoms are directly rejected. Default: True.
        n_jobs=4                        # Accelerate using 4 CPU cores. Default: 1 
    )

    # Configure deduplication and unit standardization
    conversion_rules = {('g/L', 'mg/L'): 'x * 1000', 
                        ('ug/L', 'mg/L'): 'x / 1000',}
    DP.config_deduplicator(condition_cols, data_type, method, custom_method, priority, standard_unit, conversion_rules, log_transform, dropna_conditions)
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

## Advanced Configuration

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
  - `streamlit>=1.0.0` (Required for GUI)
- Optional Dependencies (install as needed, if not installed, then send the request using `requests`.):
  - `pubchempy>=1.0.5`: For PubChem integration
  - `chemspipy>=2.0.0`: For ChemSpider (requires API key)
  - `ctx-python>=0.0.1a10`: For CompTox Dashboard (requires API key)

## License
Apache License 2.0 - See [LICENSE](LICENSE) for details

## Support
Report issues on [GitHub Issues](https://github.com/Hya0FAD/DiPTox/issues)
