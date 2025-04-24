# MolPPC - A Python Package for Molecular Preprocessing Pipeline for Chemistry

A toolkit for molecular data preprocessing, standardization, and multi-source molecular identifier integration.

## Features
- **Molecular Preprocessing**: 
  - Neutralize charges
  - Remove salts and solvents
  - Validate atomic composition
  - Remove stereochemistry
  - Hydrogen removal
  - Chemical sanitization
  - Fragment handling (keep largest fragment)
  
- **Deduplication**:
  - Flexible deduplication strategies
  - Customizable matching conditions
  
- **Chemical Identifier Integration**:
  - CAS number validation
  - IUPAC name standardization
  - SMILES processing
  
- **Web Data Integration**:
  - Multi-source chemical property queries
  - Concurrent requests support
  - API key management

## Installation
```bash
pip install -i https://test.pypi.org/simple/ molppc
```
## Quick Start
```{python}
from molppc import MolecularProcessor

# Initialize processor
MP = MolecularProcessor()

# Load data
MP.load_data(input_data='file_path/list/dataframe', smiles_col, target_col, cas_col)

# Configure preprocessing
MP.preprocess(remove_salts, remove_solvents, remove_mixtures, remove_inorganic, neutralize, check_valid_atoms, remove_stereo, remove_hs, keep_largest_fragment, hac_threshold)
# remove_salts: Whether to remove salts.
# remove_solvents: Whether to remove solvent molecules.
# remove_mixtures: Whether to remove mixtures.
# remove_inorganic: Whether to remove inorganic molecules.
# neutralize: Whether to neutralize charges.
# check_valid_atoms: Whether to check for valid atoms.
# remove_stereo: Whether to remove stereochemistry.
# remove_hs: Whether to remove hydrogen atoms.
# keep_largest_fragment: Whether to keep the largest fragment.
# hac_threshold: Threshold for salt removal (heavy atoms count).

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
  - `ctx-python>=0.0.1a7`: For CompTox Dashboard (requires API key and Python>=3.10)

## License
MIT License - See [LICENSE](LICENSE) for details

## Support
Report issues on [GitHub Issues](https://github.com/Hya0FAD/MolPPC/issues)
