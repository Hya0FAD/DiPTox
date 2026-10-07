# ==============================================================================
# DiPTox - Usage Examples & Best Practices
# ==============================================================================
# Python >=3.8. Install this source checkout from its project directory:
# python -m pip install .
# For the published release: python -m pip install --upgrade diptox
# GUI: diptox-gui (or python -m diptox.cli)
# CLI: python -m diptox --help
# Python 3.8/3.9 installs NiceGUI 2.24.2; Python >=3.10 installs NiceGUI >=3.16,<4.
# Both use the same current GUI and processing functions.
# Run this file for the offline SMILES example. Other examples need your input files.

# ⚠️ CRITICAL NOTE FOR WINDOWS USERS:
# When using multi-processing (n_jobs > 1) in 'preprocess', you MUST wrap your
# code logic inside an `if __name__ == '__main__':` block.
# Failure to do so will result in recursive process loops and memory explosion.


import multiprocessing

import pandas as pd

from diptox import DiptoxPipeline

# ==============================================================================
# Use Case 1: Full Workflow for File Processing
# ==============================================================================
# This example shows a complete pipeline from loading a data file to saving the processed results.
def use_case_1(input_path="path/to/your/FileName.xlsx",
               output_path="Processed_FileName.csv", enrich=False):
    a = DiptoxPipeline(interactive=False)
    a.load_data(input_data=str(input_path), # Supports .xls, .csv, .txt, .sdf, .mol, .smi
                smiles_col='SMILES_Column_Name', 
                target_col='Target_Column_Name', 
                unit_col='Unit_Column_Name',
                cas_col='CAS_Column_Name', 
                id_col='ID_Column_Name',
                header=0) # Map your actual columns; omit roles your dataset does not have.
    # For multi-sheet Excel files, also pass sheet_name='YourSheet' to load_data.
    a.preprocess(mixture_mode='reject', element_policy='allow_all', n_jobs=1)
    # Choose the MW basis matching the measured substance. For reported MW, use
    # molecular_weight_col='Reported MW' instead of molecular_weight_source.
    a.standardize_units(standard_unit='mg/L', molecular_weight_source='original')
    a.config_deduplicator(data_type='continuous', # Numeric endpoint values in this example.
                        method='IQR', aggregation='mean',
                        condition_cols=['pH', 'temperature'], # These columns must exist.
                        dropna_conditions=False,
                        log_transform='None')
    # Categorical endpoints: skip unit conversion and use data_type='discrete',
    # method='vote' (or method='priority', priority=['Active', 'Inactive']).
    # Use '-log10' only for an appropriate positive numeric endpoint/scale.
    a.dataset_deduplicate()
    if enrich: # Optional: sends SMILES to PubChem; requires network access.
        a.config_web_request(sources=['pubchem'], max_workers=2, retries=1)
        a.web_request(send='smiles', request=['cas', 'iupac'])
    a.save_results(str(output_path)) # Supports .xlsx, .csv, .txt, .sdf, .smi
    print(a.get_history())
    return a


# ==============================================================================
# Use Case 2: Full Workflow for a List of SMILES
# ==============================================================================
# This example demonstrates how to process a simple list of SMILES strings and customize the chemical processing rules.
def use_case_2(output_path='Processed_SMILES.csv', n_jobs=1, enrich=False):
    b = DiptoxPipeline(interactive=False)
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
        r'[H][C@]12O/C=C\[C@@]1([H])C5=C(O2)C=C(OC)C4=C5OC=3C=CC=C(O)C=3C4=O'
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
        mixture_mode='largest',    # 'keep', 'reject' (default), or 'largest'.
        hac_threshold=3,            # Minimum heavy atoms for largest-fragment selection.
        remove_inorganic=False,     # Remove common inorganic molecules. Default: True.
        neutralize=True,            # Neutralize charges on the molecule. Default: True.
        reject_non_neutral=False,   # Only retain the molecules whose formal charge is zero. Default: False.
        element_policy='allowed_atoms', # Reject whole molecules with unlisted atoms.
                                   # Alternatives: 'allow_all' (default), 'reject_metals'.
        remove_stereo=False,        # Remove stereochemistry information. Default: False.
        remove_isotopes=True,       # Remove isotopic information. Default: True.
        remove_hs=True,             # Remove explicit hydrogen atoms. Default: True.
        add_hs=False,               # If True, add explicit H after other processing.
        sanitize=True,
        reject_radical_species=True, # Reject radical species. Default: True.
        n_jobs=n_jobs,              # Use >1 inside the guarded main block on Windows.
    )

    b.filter_by_atom_count(min_total_atoms=4)
    # min_heavy_atoms
    # max_heavy_atoms
    # min_total_atoms
    # max_total_atoms
    b.config_deduplicator(data_type='smiles') # Explicit structure-only deduplication.
    b.dataset_deduplicate()
    b.substructure_search(query_pattern='[C@@H]', is_smarts=True)  # Find substructures
    if enrich: # Optional: sends SMILES to PubChem.
        b.config_web_request(sources=['pubchem'], max_workers=2, retries=1)
        b.web_request(send='smiles', request=['cas', 'iupac'])
    b.save_results(str(output_path))
    print(b.df['Canonical SMILES'])
    print(b.get_history())
    # Rejected/filtered records remain available for auditing:
    print(b.get_excluded_records())
    # b.save_excluded_results('Excluded_SMILES.csv')
    # b.undo()  # Restore the previous pipeline checkpoint, if needed.
    return b


# ==============================================================================
# Use Case 3: Fetch SMILES from CAS, then Process
# ==============================================================================
# This example shows how to use CAS numbers as the primary input, fetch SMILES from the web, and then process the data.
def use_case_3(input_path="path/to/your/FileName.csv",
               output_path='CAS_Processed_Results.xlsx'):
    # This example sends CAS identifiers to PubChem and requires network access.
    c = DiptoxPipeline(interactive=False)
    c.load_data(input_data=str(input_path),
                cas_col='CAS_Number_Column', target_col='Value_Column', unit_col='Unit_Column')
    c.config_web_request(sources=['pubchem'], max_workers=2, retries=1)
    c.web_request(send='cas', request=['smiles', 'iupac'])
    c.preprocess()
    c.substructure_search(query_pattern='C(=O)O', is_smarts=True)  # Match fragments
    rules = {
        ('ng/mL', 'mg/L'): 'x / 1000',
        ('10^-6 M', 'mg/L'): 'x * mw / 1000'
    }
    c.config_deduplicator(data_type='continuous',        # Defaults to 'auto' method
                        standard_unit='mg/L',          # Target unit
                        conversion_rules=rules,        # Direct rules to the target unit.
                        molecular_weight_source='original',
                        aggregation='mean',
                        log_transform='None')
    c.dataset_deduplicate()
    c.save_results(str(output_path))
    return c


# ==============================================================================
# Use Case 4: SMILES-only Deduplication (for Pre-training)
# ==============================================================================
# A simple workflow for when you only need to standardize and deduplicate a large list of SMILES.
def use_case_4(input_path="path/to/your/Large_SMILES_Dataset.csv",
               output_path='Deduplicated_SMILES.csv', n_jobs=1):
    d = DiptoxPipeline(interactive=False)
    d.load_data(input_data=str(input_path),
                smiles_col='Smiles')
    d.preprocess(mixture_mode='reject', n_jobs=n_jobs)
    d.config_deduplicator(data_type='smiles')
    d.dataset_deduplicate()
    d.save_results(str(output_path))
    print(f"Use Case 4 finished: {output_path}")
    return d


# ==============================================================================
# Use Case 5: Handling SDF/SMI Files
# ==============================================================================
# This example demonstrates loading data from .sdf or .smi files and saving in various formats.
def use_case_5(input_path="path/to/your/zinc_database.smi",
               output_path='Processed_Molecules.sdf'):
    e = DiptoxPipeline(interactive=False)

    # Example for loading from an SDF file:
    # e.load_data(input_data=r"path/to/your/AllPublicnew.sdf",
    #             smiles_col='SMILES', target_col='ReadyBiodegradability', cas_col='CASRN')

    # Example for loading from an SMI file:
    # A headerless two-column SMI file: SMILES followed by a whitespace-separated ID.
    e.load_data(input_data=str(input_path),
                smiles_col='smiles', id_col='zinc_id', header=None)
                
    e.preprocess()
    e.config_deduplicator(data_type='smiles') # No categorical target in this file.
    e.dataset_deduplicate()
    e.save_results(str(output_path))
    print(f"Use Case 5 finished: {output_path}")
    return e


# ==============================================================================
# Advanced Use Case: Custom Outlier Detection Logic
# ==============================================================================
# You can provide your own function to handle outlier detection during continuous data deduplication.
def custom_outlier_filter(values: pd.Series):
    """
    A custom filter that keeps values in the range [-3.5, 3.5].
    Returns a filtered Series with its original index and a method description.
    The deduplicator applies the configured aggregation afterward.
    """
    clean_values = values[(values >= -3.5) & (values <= 3.5)]
    
    if len(clean_values) >= 1:
        return clean_values, "custom_range"
    else: 
        # Explicit example policy: retain the original group if all are outside.
        return values, "fallback_keep_all"


def advanced_use_case(input_path="path/to/your/FileName.xlsx",
                      output_path='Custom_Deduplication_Results.csv'):
    f = DiptoxPipeline(interactive=False)
    f.load_data(input_data=str(input_path),
                smiles_col='Smiles', target_col='Value')
                
    f.preprocess(element_policy='allow_all')
    f.config_deduplicator(custom_method=custom_outlier_filter, data_type='continuous',
                         aggregation='mean', log_transform='None')
    f.dataset_deduplicate()

    f.save_results(str(output_path))
    return f


# ==============================================================================
# Use Case 6: Filter Rows by Column Values
# ==============================================================================
def use_case_6():
    """Discover column values and filter rows; omitted selections keep all rows."""
    pipeline = DiptoxPipeline(interactive=False)
    pipeline.load_data(pd.DataFrame({
        'SMILES': ['CCO', 'CCC', 'CCCC', 'CC'],
        'Species': ['rat', 'mouse', 'rat', None],
        'Study type': ['in vivo', 'in vitro', 'in vitro', 'in vivo'],
    }), smiles_col='SMILES')
    print(pipeline.get_column_values('Species'))  # [{'value': 'rat', 'count': 2}, ...]
    pipeline.filter_by_values('Species', ['rat'], mode='keep')
    pipeline.filter_by_values('Study type', ['in vitro'], mode='remove')
    # None inside a selection matches missing cells: values=[None].
    # values=[] (or omitted values) does not remove any rows, in either mode.
    print(pipeline.df)
    print(pipeline.get_excluded_records())
    # pipeline.save_excluded_results('Column_filter_excluded.csv')
    # pipeline.undo()  # Undo the latest column filter, including its exclusions.
    return pipeline


def use_case_7():
    """Convert exposure duration before using it as a deduplication condition."""
    pipeline = DiptoxPipeline(interactive=False)
    pipeline.load_data(pd.DataFrame({
        'SMILES': ['CCO', 'CCO'], 'Value': [1., 3.], 'Unit': ['mg/L', 'mg/L'],
        'Duration': [1., 24.], 'Duration unit': ['d', 'h'],
    }), smiles_col='SMILES', target_col='Value', unit_col='Unit')
    value_col, unit_col = pipeline.transform_column(
        value_col='Duration', unit_col='Duration unit', standard_unit='h',
        log_transform='None',  # 'log10' or '-log10'; conversion runs first.
    )
    pipeline.config_deduplicator(condition_cols=[value_col, unit_col],
                                 data_type='continuous', aggregation='mean')
    pipeline.dataset_deduplicate()
    print(pipeline.df)
    return pipeline


def use_case_8():
    """Combine categories into one condition, retaining the original column."""
    pipeline = DiptoxPipeline(interactive=False)
    pipeline.load_data(pd.DataFrame({
        'SMILES': ['CCO', 'CCO', 'CCO'],
        'Species': ['mouse', 'rabbit', 'human'], 'Value': [1., 3., 5.],
    }), smiles_col='SMILES', target_col='Value')
    print(pipeline.get_column_values('Species'))
    group = pipeline.merge_column_values(
        'Species', ['mouse', 'rabbit'], replacement='other', output_column='Species group')
    # The new column is ['other', 'other', 'human']; source values and rows remain intact.
    pipeline.config_deduplicator(condition_cols=[group], data_type='continuous', aggregation='mean')
    pipeline.dataset_deduplicate()
    print(pipeline.df)
    return pipeline


def use_case_9():
    """Apply several category mappings simultaneously into one output column."""
    import json
    from pathlib import Path
    rules = json.loads((Path(__file__).parent / 'examples' / 'life_stage_merge_groups.json').read_text(encoding='utf-8'))
    pipeline = DiptoxPipeline(interactive=False)
    pipeline.load_data(pd.DataFrame({'SMILES': ['CCO'] * 6,
                                    'Stage': ['Egg', 'Fry', 'Juvenile', 'Mature', 'Sperm', 'Unknown']}),
                       smiles_col='SMILES')
    output = pipeline.merge_column_values('Stage', groups=rules, output_column='Stage group')
    print(pipeline.df[['Stage', output]])
    # Select output as a condition column when configuring deduplication.
    return pipeline


# Main Execution Block (required for Windows multiprocessing)
if __name__ == '__main__':
    multiprocessing.freeze_support()
    # Uncomment the function you want to run:
    
    # use_case_1()
    use_case_2()
    # use_case_3()
    # use_case_4()
    # use_case_5()
    # use_case_6()
    # use_case_7()
    # use_case_8()
    # use_case_9()
    # advanced_use_case()
