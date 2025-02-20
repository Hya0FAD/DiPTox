# MolPPC - A Python Package for Molecular Data Preprocessing
A toolkit for batch querying chemical properties from multiple web sources.

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
  - SMILES/SMARTS processing
  
- **Web Data Integration**:
  - Multi-source chemical property queries
  - Concurrent requests support
  - API key management

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
MP.process(neutralize, remove_salts, check_valid_atoms, remove_stereo, remove_hs, keep_largest_fragment, hac_threshold, sanitize)

# Configure deduplication
MP.config_deduplicator(condition_cols, data_type, method, custom_method)
MP.deduplicate_data()

# Configure web queries
MP.config_web_request(source='pubchem/chemspider/comptox/cactus', max_workers, ...)
MP.add_web_request(cas, iupac, smiles)

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
- Optional Dependencies (install as needed):
  - `pubchempy>=1.0.4`: For PubChem integration
  - `chemspipy>=2.0.0`: For ChemSpider (requires API key)
  - `ctx-python>=0.0.1a7`: For CompTox Dashboard (requires API key)

## License
MIT License - See [LICENSE](LICENSE) for details

## Support
Report issues on [GitHub Issues](https://github.com/Hya0FAD/MolPPC/issues)
