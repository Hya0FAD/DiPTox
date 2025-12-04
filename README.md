# DiPTox - Data Integration and Processing for Computational Toxicology

![PyPI Test Version](https://img.shields.io/badge/testpypi-1.3.1-blue) ![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg) ![Python Version](https://img.shields.io/badge/python-3.8+-brightgreen.svg) [![Chinese](https://img.shields.io/badge/-%E4%B8%AD%E6%96%87%E7%89%88-blue.svg)](./README_ZH.md)
<p align="center">
  <img src="assets/TOC.png" alt="DiPTox Workflow Diagram" width="500">
</p>
**DiPTox** is a Python toolkit designed for the robust preprocessing, standardization, and multi-source data integration of molecular datasets, with a focus on computational toxicology workflows.

## New in v1.3: Unit Standardization & Transformation
Handling heterogeneous experimental data often involves dealing with messy units (e.g., *mg/L, ug/mL, M, ppm, %*). DiPTox v1.3 introduces **Unit Processor** to automate this process:
-   **Automatic Conversion**: Built-in rules for **Concentration** (mass/vol, molar, parts-per), **Time**, **Pressure**, and **Temperature**.
-   **Custom Formulas**: Define your own mathematical rules (e.g., `x * 1000` or `10**(-x)`) interactively via the GUI or script.
-   **Log Transformation**: The deduplication module now supports optional `-log10` transformation (e.g., converting IC50 to pIC50) with a single parameter.
-   **Enhanced Inorganic Filtering**: Significantly improved the `remove_inorganic` module using strict SMARTS pattern matching. It now accurately correctly identifies and removes complex inorganic species (e.g., ionic cyanides `[C-]#N`, carbonates, carbonyls) without misclassifying organic structures (like nitriles).
-   **Step-by-Step Audit Log**: Introduced a comprehensive **History Tracking** system. DiPTox now automatically records a timeline of every operation (Loading, Preprocessing, Filtering, Deduplication, etc.), tracking the **timestamp**, **operation details**, and the **number of rows retained or removed (Delta)** at each stage. This feature is available in both the Python API (`get_history()`) and the GUI.

## DiPTox Community Check-in (Optional)
To help us understand our user base and improve the software, DiPTox includes a one-time, optional survey on first use. 
-   **Completely Optional**: You can skip it with a single click.
-   **Privacy-Focused**: The information helps us with academic impact assessment. It will not be shared.

## Core Features

#### Graphical User Interface (GUI)
Powered by Streamlit, the GUI allows users to perform all workflows visually without writing code.
-   **Visual Operation**: Complete workflow control via a web browser.
-   **Real-time Preview**: Instantly view data changes after applying rules.
-   **Rule Management**: Add/Remove valid atoms, salts, solvents, and **unit conversion formulas** interactively.

#### Chemical Preprocessing & Standardization
A configurable pipeline to clean and normalize chemical structures in a specific, controlled order:
-   **Remove salts**
-   **Remove solvents**
-   **Handle mixtures** (e.g., keep the largest fragment)
-   **Remove inorganic molecules**
-   **Neutralize charges**
-   **Validate atomic composition** against a list of allowed elements
-   **Remove explicit hydrogens**
-   **Remove stereochemistry**
-   **Remove isotopes**
-   **Handle Radical**
-   **Standardize molecules** to canonical SMILES
-   **Filter by atom count** (heavy or total atoms)

#### Unit Standardization (New)
-   Standardize diverse target values to a single unit (e.g., converting all data to `mg/L`).
-   Support for complex conversions and custom math expressions.

#### Data Deduplication
-   Flexible strategies for handling duplicate entries (`smiles` or `continuous`/`discrete` targets).
-   Customizable matching conditions (e.g., temperature, pressure) and deduplication methods (`auto`, `IQR`, `3sigma`, or custom functions).
-   **New:** Optional `-log10` transformation for bioactivity data before deduplication.

#### Identifier & Property Integration via Web Services
-   Fetch and interconvert chemical identifiers (**CAS, SMILES, IUPAC Name, Common Name, MW**) from multiple online databases (**PubChem, ChemSpider, CompTox, Cactus, CAS Common Chemistry, ChEMBL**).
-   High-performance **concurrent requests** to accelerate data retrieval.
-   Centralized API key management for services that require authentication.

#### Utility Tools
-   Perform **substructure searches** using SMILES or SMARTS patterns.
-   **Customize chemical processing rules** for neutralization reactions, salt/solvent lists, and valid atoms.
-   **Display a summary** of all currently active processing rules.

## Installation
```bash
pip install -i https://test.pypi.org/simple/ diptox
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
  remove_salts=True,          # Remove salt fragments. Default: True.
  remove_solvents=True,       # Remove solvent fragments. Default: True.
  remove_mixtures=True,       # Handle mixtures based on fragment size. Default: False.
  hac_threshold=3,            # Heavy atom count threshold for fragment removal. Default: 3.
  keep_largest_fragment=True, # Keep the largest fragment in a mixture. Default: True.
  remove_inorganic=False,     # Remove common inorganic molecules. Default: True.
  neutralize=True,            # Neutralize charges on the molecule. Default: True.
  reject_non_neutral=False,   # Only retain the molecules whose formal charge is zero. Default: False.
  check_valid_atoms=True,     # Check if all atoms are in the valid list. Default: False.
  strict_atom_check=False,    # If True, discard molecules with invalid atoms. If False, try to remove them. Default: False.
  remove_stereo=False,        # Remove stereochemistry information. Default: False.
  remove_isotopes=True,       # Remove isotopic information. Default: True.
  remove_hs=True,             # Remove explicit hydrogen atoms. Default: True.
  reject_radical_species=True # Molecules containing free radical atoms are directly rejected. Default: True.
)

# Configure deduplication and unit standardization
conversion_rules = {('g/L', 'mg/L'): 'x * 1000', 
                    ('ug/L', 'mg/L'): 'x / 1000',}
DP.config_deduplicator(condition_cols, data_type, method, custom_method, priority, standard_unit, conversion_rules, log_transform)
DP.data_deduplicate()

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
