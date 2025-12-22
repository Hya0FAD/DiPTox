# ==============================================================================
# Installation & Requirements
# ==============================================================================

# To install from TestPyPI:
# pip install -i https://test.pypi.org/simple/ diptox

# Requirements:
# - python>=3.8
# - rdkit>=2023.3
# - requests
# - tqdm
# - scipy
# - openpyxl

# Optional Dependencies:
# - pubchempy (for PubChem integration)
# - ctx-python (for CompTox integration)
# - chemspipy (for ChemSpider integration)

# Latest Version: 1.3.5 (2025.12.22)


# ==============================================================================
# Use Case 1: Full Workflow for File Processing
# ==============================================================================
# This example shows a complete pipeline from loading a data file to saving the processed results.

from diptox import DiptoxPipeline

a = DiptoxPipeline()
a.load_data(input_data=r"path/to/your/FileName.xlsx", # Supports .xls, .csv, .txt, .sdf, .smi
            smiles_col='SMILES_Column_Name', 
            target_col='Target_Column_Name', 
            unit_col='Unit_Column_Name',
            cas_col='CAS_Column_Name', 
            id_col='ID_Column_Name',
            header=0) # Note: smiles_col, target_col, cas_col, id_col are all optional.
a.preprocess()
a.standardize_units(standard_unit='mg/L')
a.config_deduplicator(data_type='discrete', # 'smiles', 'discrete' or 'continuous'
                      method='vote', # 'auto', 'vote', '3sigma', 'IQR'
                      condition_cols=['pH','temperature'],
                      log_transform=True) # Example: provide a list of condition columns
a.data_deduplicate()
a.config_web_request(sources=['pubchem', 'cas', 'comptox', 'chemspider', 'cactus'],
                     comptox_api_key='your_key_here', 
                     chemspider_api_key='your_key_here', 
                     cas_api_key='your_key_here', 
                     max_workers=4)
a.web_request(send='smiles', request=['cas', 'iupac'])
a.save_results(r"path/to/your/Processed_FileName.csv") # Can save as .xlsx, .csv, .txt, .sdf, .smi


# ==============================================================================
# Use Case 2: Full Workflow for a List of SMILES
# ==============================================================================
# This example demonstrates how to process a simple list of SMILES strings and customize the chemical processing rules.

from diptox import DiptoxPipeline

b = DiptoxPipeline()
smiles = [
    "C(C(=O)O)C.C(C(=O)O)C.[Na+]", 
    "C1=CC=CC=C1.CCCC.[Hg+2].[Na+]",
    "CC1=CC=CC=C1.CC1=CC=CC=C1.CC(O)C", 
    "[Hg+2].[Cl-].[Na+]", 
    "Cl[Zr](Cl)=O",
    'CCC(=O)O.CCO', 
    "C=O", 
    'CCCC', 
    'CC', 
    'CCC.CCCC', 
    '234', 
    '', 
    'C(Cl)(Cl)(Cl)Cl',
    'NaCl', 
    '[Si](C)(C)[Si](C)(C)C', 
    '[N-]C(=O)C',
    'OC[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)CO',
    'C1=CC(=C2C(=C1)OC(O2)(F)F)C3=CNC=C3C#N', 
    'C#C', 
    'OC(C(O)C(=O)O)C(=O)O.CCC',
    'CO', 
    'O=C1CN(/N=C/c2ccc([N+](=O)[O-])o2)C(=O)N1',
    '[H][C@]12O/C=C\[C@@]1([H])C5=C(O2)C=C(OC)C4=C5OC=3C=CC=C(O)C=3C4=O'
]
b.load_data(input_data=smiles, smiles_col='SMILES')

# Customize chemical processing rules
b.remove_neutralization_rule('[$([N-]C=O)]')          # Remove a default neutralization rule
b.add_neutralization_rule('[$([N-]C=O)]', 'N')        # Add a new neutralization rule
b.manage_atom_rules(atoms=['Si', 'Zr'], add=True)     # Add/remove atoms for validation
b.manage_default_salt(salts=['[Hg+2]', '[Ba+2]'], add=True) # Add/remove salts from the salt list
b.manage_default_solvent(solvents='CCC', add=True)    # Add a custom solvent

b.display_processing_rules()

b.preprocess(
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

b.filter_by_atom_count(min_total_atoms=4)
# min_heavy_atoms
# max_heavy_atoms
# min_total_atoms
# max_total_atoms
b.config_deduplicator()  # With no arguments, defaults to SMILES-only deduplication
b.data_deduplicate()
b.substructure_search(query_pattern='[C@@H]', is_smarts=True)  # Find substructures
b.config_web_request(sources=['pubchem'], max_workers=4)
b.web_request(send='smiles', request=['cas', 'iupac'])
b.save_results('Processed_SMILES.csv')
print(b.df['Canonical smiles'])


# ==============================================================================
# Use Case 3: Fetch SMILES from CAS, then Process
# ==============================================================================
# This example shows how to use CAS numbers as the primary input, fetch SMILES from the web, and then process the data.

from diptox import DiptoxPipeline

c = DiptoxPipeline()
c.load_data(input_data="path/to/your/FileName.csv",
            cas_col='CAS_Number_Column', target_col='Value_Column')
c.config_web_request(sources=['comptox'], comptox_api_key='your_key_here', max_workers=4)
c.web_request(send='cas', request=['smiles', 'iupac'])
c.preprocess()
c.substructure_search(query_pattern='C(=O)O', is_smarts=True)  # Match fragments
rules = {
    ('ng/mL', 'mg/L'): 'x / 1000000',
    ('10^-6 M', 'M'): 'x * 1e-6'
}
c.config_deduplicator(data_type='continuous',        # Defaults to 'auto' method
                      standard_unit='mg/L',          # Target unit
                      conversion_rules=rules,      # Custom rules
                      log_transform=False)  
c.data_deduplicate()
c.save_results('CAS_Processed_Results.xlsx')


# ==============================================================================
# Use Case 4: SMILES-only Deduplication (for Pre-training)
# ==============================================================================
# A simple workflow for when you only need to standardize and deduplicate a large list of SMILES.

from diptox import DiptoxPipeline

d = DiptoxPipeline()
d.load_data(input_data="path/to/your/Large_SMILES_Dataset.csv",
            smiles_col='Smiles')
d.preprocess()
d.config_deduplicator()
d.data_deduplicate(data_type='smiles')
d.save_results('Deduplicated_SMILES.csv')
print("Use Case 4 finished. Check for 'Deduplicated_SMILES.csv'.")


# ==============================================================================
# Use Case 5: Handling SDF/SMI Files
# ==============================================================================
# This example demonstrates loading data from .sdf or .smi files and saving in various formats.

from diptox import DiptoxPipeline

e = DiptoxPipeline()

# Example for loading from an SDF file:
# e.load_data(input_data=r"path/to/your/AllPublicnew.sdf",
#             smiles_col='SMILES', target_col='ReadyBiodegradability', cas_col='CASRN')

# Example for loading from an SMI file:
e.load_data(input_data=r"path/to/your/zinc_database.smi",
            smiles_col='smiles', id_col='zinc_id')
            
e.preprocess()
e.config_deduplicator(data_type='discrete')
e.data_deduplicate()
e.config_web_request(sources=['pubchem'], max_workers=4)
e.web_request(send='smiles', request=['cas'])

e.save_results(r"Processed_Molecules.sdf")
print("Use Case 5 finished. Check for 'Processed_Molecules.sdf'.")


# ==============================================================================
# Advanced Use Case: Custom Outlier Detection Logic
# ==============================================================================
# You can provide your own function to handle outlier detection during continuous data deduplication.

from diptox import DiptoxPipeline
import pandas as pd
import numpy as np

def custom_outlier_filter(values: pd.Series):
    """
    A custom filter that keeps values in the range [-3.5, 3.5].
    Returns the mean of the filtered values and a description of the method.
    """
    clean_values = values[(values >= -3.5) & (values <= 3.5)]
    
    if len(clean_values) >= 1:
        # If there are any clean values, return their mean
        return clean_values.mean(), "custom_range_mean"
    else: 
        # If all values are outliers, fall back to the mean of original values
        return values.mean(), "fallback_mean"

f = DiptoxPipeline()
f.load_data(input_data="path/to/your/FileName.xlsx",
            smiles_col='Smiles', target_col='Value')
            
f.preprocess(check_valid_atoms=False)
f.config_deduplicator(custom_method=custom_outlier_filter, data_type='continuous')
f.data_deduplicate()

f.save_results('Custom_Deduplication_Results.csv')