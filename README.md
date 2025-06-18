# MolPPC - A Python Package for Molecular Preprocessing Pipeline for Chemistry

A toolkit for molecular data preprocessing, standardization, and multi-source molecular identifier integration.

## Features
- **Molecular Preprocessing**: 
  - Neutralize charges
  - Remove salts and solvents
  - Validate atomic composition
  - Filter by atom count (heavy or total)
  - Remove stereochemistry
  - Hydrogen removal
  - Chemical sanitization
  - Fragment handling (keep largest fragment)
  
- **Deduplication**:
  - Flexible deduplication strategies (`smiles_only`, `continuous`, `discrete`)
  - Customizable matching conditions and outlier handling methods (`auto`, `IQR`, `3sigma`)
  
- **Chemical Identifier Integration**:
  - CAS number validation
  - IUPAC name standardization
  - SMILES processing
  
- **Web Data Integration**:
  - Multi-source chemical property queries (PubChem, ChemSpider, CompTox, Cactus)
  - Concurrent requests support
  - API key management

## What's New in Version 0.12.0

This release introduces several new features and significant improvements to enhance flexibility and user control over the preprocessing pipeline.

### ✨ New Features

1.  **Atom Count Filtering**:
    - A new method `filter_by_atom_count()` has been added to `MolecularProcessor`.
    - Users can now filter the dataset based on minimum or maximum counts of **heavy atoms** or **total atoms** (including hydrogens). This is useful for removing fragments, small molecules, or overly large structures.
    - Example:
      ```{python}
      # Keep molecules with 5 to 50 heavy atoms
      MP.filter_by_atom_count(min_heavy_atoms=5, max_heavy_atoms=50)

      # Keep molecules with a total atom count (including H) of at least 10
      MP.filter_by_atom_count(min_total_atoms=10)
      ```

2.  **Enhanced Chemical Rules Management**:
    - Added `display_processing_rules()` to `MolecularProcessor`, which provides a clear, formatted summary of all currently active chemical rules.
    - This includes valid atoms, neutralization rules, and the dynamically generated lists of effective salts and solvents (respecting all user additions and removals).
    - Example:
      ```{python}
      MP.manage_default_solvent("c1ccccc1", add=False) # Remove benzene
      MP.display_processing_rules() # See the updated list of rules
      ```

### ♻️ API Improvements & Refinements

1.  **Simplified Deduplication `data_type`**:
    - The `DataDeduplicator` has been refactored for clarity. To perform simple deduplication based only on SMILES (and optional condition columns), you now set `data_type='smiles'`.
    - Example:
      ```{python}
      # New, clearer way to perform simple SMILES-based deduplication
      MP.config_deduplicator(data_type='smiles')
      MP.data_deduplicate()
      ```

2.  **Robust Salt & Solvent Handling**:
    - The internal logic for managing salt and solvent lists has been completely overhauled.
    - Custom-added salts/solvents that already exist in the default lists are now correctly handled to avoid duplication.
    - The `SaltRemover` instance from RDKit is now guaranteed to use **only** the lists defined by `MolPPC` and user modifications, completely ignoring RDKit's built-in defaults. This gives users full control.

### 🐛 Bug Fixes

-   Fixed a critical bug in atom count calculation where `GetNumAtoms()` would incorrectly report only heavy atoms if hydrogens were not explicit. The `filter_by_atom_count` feature now correctly calculates total atoms.

---

## Installation
```bash
pip install molppc
```
## Quick Start
```{python}
from molppc import MolecularProcessor

# Initialize processor
MP = MolecularProcessor()

# Load data
MP.load_data(input_data='file_path/list/dataframe', smiles_col, target_col, cas_col)

# Configure preprocessing
MP.preprocess(neutralize, remove_salts, check_valid_atoms, remove_stereo, remove_hs, keep_largest_fragment, hac_threshold, sanitize)

# Configure deduplication
MP.config_deduplicator(condition_cols, data_type, method, custom_method)
MP.data_deduplicate()

# Configure web queries
MP.config_web_request(source='pubchem/chemspider/comptox/cactus', max_workers, ...)
MP.web_request(cas, iupac, smiles)

# Substructure search
MP.substructure_search(query_pattern, is_smarts=True)

# Save results
MP.save_results(output_path='file_path')
```
## Advanced Configuration
### Web Service Integration
- `PubChem`: https://pubchem.ncbi.nlm.nih.gov/
- `ChemSpider`: https://www.chemspider.com/
- `CompTox`: https://comptox.epa.gov/dashboard/
- `Cactus`: https://cactus.nci.nih.gov/
### Note: Add API keys during initialization:
```{python}
MP.config_web_request(
    chemspider_api_key='your_key',
    comptox_api_key='your_key'
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
Report issues on [GitHub Issues](https://github.com/Hya0FAD/MolPPC/issues)
