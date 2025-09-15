# DiPTox - Data Integration and Processing for Computational Toxicology

![PyPI Test Version](https://img.shields.io/badge/testpypi-1.0.1-blue) ![License](https://img.shields.io/badge/license-MIT-blue.svg) ![Python Version](https://img.shields.io/badge/python-3.8+-brightgreen.svg) [![Chinese](https://img.shields.io/badge/-%E4%B8%AD%E6%96%87%E7%89%88-blue.svg)](./README_ZH.md)
<p align="center">
  <img src="assets/TOC.png" alt="DiPTox Workflow Diagram" width="500">
</p>
**DiPTox** is a Python toolkit designed for the robust preprocessing, standardization, and multi-source data integration of molecular datasets, with a focus on computational toxicology workflows.

## Core Features

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

#### Data Deduplication
-   Flexible strategies for handling duplicate entries (`smiles` or `continuous`/`discrete` targets).
-   Customizable matching conditions (e.g., temperature, pressure) and deduplication methods (`auto`, `IQR`, `3sigma`, or custom functions).

#### Identifier & Property Integration via Web Services
-   Fetch and interconvert chemical identifiers (**CAS, SMILES, IUPAC Name**) from multiple online databases (**PubChem, ChemSpider, CompTox, Cactus**).
-   High-performance **concurrent requests** to accelerate data retrieval.
-   Centralized API key management for services that require authentication.

#### Utility Tools
-   Perform **substructure searches** using SMILES or SMARTS patterns.
-   **Customize chemical processing rules** for neutralization reactions, salt/solvent lists, and valid atoms.
-   **Display a summary** of all currently active processing rules.

## Installation
```bash
pip install diptox
```

## Quick Start
```{python}
from diptox import DiptoxPipeline

# Initialize processor
DP = DiptoxPipeline()

# Load data
DP.load_data(input_data='file_path/list/dataframe', smiles_col, target_col, cas_col)

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

# Configure deduplication
DP.config_deduplicator(condition_cols, data_type, method, custom_method)
DP.data_deduplicate()

# Configure web queries
DP.config_web_request(source='pubchem/chemspider/comptox/cactus', max_workers, ...)
DP.web_request(send='cas', request=['smiles', 'iupac'])

# Substructure search
DP.substructure_search(query_pattern, is_smarts=True)

# Save results
DP.save_results(output_path='file_path')
```

## Advanced Configuration

### Web Service Integration
DiPTox supports the following chemical databases:
-   `PubChem`: https://pubchem.ncbi.nlm.nih.gov/
-   `ChemSpider`: https://www.chemspider.com/
-   `CompTox`: https://comptox.epa.gov/dashboard/
-   `Cactus`: https://cactus.nci.nih.gov/

**Note:** `ChemSpider` and `CompTox` require API keys. Provide them during configuration:
```{python}
DP.config_web_request(
    source='chemspider/comptox',
    chemspider_api_key='your_personal_key',
    comptox_api_key='your_personal_key'
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
- Optional Dependencies (install as needed):
  - `pubchempy>=1.0.4`: For PubChem integration
  - `chemspipy>=2.0.0`: For ChemSpider (requires API key)
  - `ctx-python>=0.0.1a7`: For CompTox Dashboard (requires API key)

## License
MIT License - See [LICENSE](LICENSE) for details

## Support
Report issues on [GitHub Issues](https://github.com/Hya0FAD/DiPTox/issues)
