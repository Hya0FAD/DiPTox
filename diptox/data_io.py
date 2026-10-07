# diptox/data_io.py
import pandas as pd
from typing import Union, List, Optional, Dict, Any
import os
import tempfile
import json
from .logger import log_manager
logger = log_manager.get_logger(__name__)


class DataHandler:
    """Data loading and saving"""

    # These names assert a source structure; derived columns such as
    # "Original Canonical SMILES" are audited only when explicitly mapped.
    SDF_SMILES_PROPERTIES = frozenset({
        'smiles', 'smile', 'smi', 'smiles_string', 'smiles string',
        'isomeric_smiles', 'isomeric smiles', 'isomericsmiles',
    })

    @staticmethod
    def _sdf_structure_key(mol) -> str:
        """Compare molecular identity without maps or ordinary explicit H."""
        from rdkit import Chem

        comparable = Chem.Mol(mol)
        for atom in comparable.GetAtoms():
            atom.SetAtomMapNum(0)
        comparable = Chem.RemoveHs(comparable)
        Chem.AssignStereochemistry(comparable, cleanIt=True, force=True)
        return Chem.MolToSmiles(comparable, canonical=True, isomericSmiles=True)

    @staticmethod
    def load_data(input_data: Union[str, List[str], pd.DataFrame],
                  smiles_col: str = None,
                  cas_col: Optional[str] = None,
                  name_col: Optional[str] = None,
                  target_col: Optional[str] = None,
                  unit_col: Optional[str] = None,
                  inchikey_col: Optional[str] = None,
                  id_col: Optional[str] = None,
                  interactive: bool = True,
                  **kwargs) -> pd.DataFrame:
        """
        Unified data loading entry point
        :param input_data: Supports three i nput types:
            - File path (csv/xlsx/xls/txt)
            - SMILES list
            - Pre-loaded DataFrame
        :param smiles_col: SMILES column name (optional)
        :param cas_col: The column name for CAS Numbers (optional)
        :param name_col: Name column name (optional)
        :param target_col: Target value column name (optional)
        :param unit_col: Unit for target value column name (optional)
        :param inchikey_col: Inchikey column name (optional)
        :param id_col: SMI file's SMILES ID column name (optional)
        :param interactive: Whether missing file choices may prompt on stdin.
        """
        if isinstance(input_data, str):
            return DataHandler._load_from_file(input_data, smiles_col, cas_col, name_col, target_col, unit_col,
                                               inchikey_col, id_col, interactive=interactive, **kwargs)
        elif isinstance(input_data, list):
            return DataHandler._load_from_list(input_data, smiles_col)
        elif isinstance(input_data, pd.DataFrame):
            return DataHandler._load_from_dataframe(input_data, smiles_col, target_col)
        else:
            logger.error(f"Unsupported input types: {type(input_data)}")

    @staticmethod
    def _load_from_file(file_path: str,
                        smiles_col: Optional[str] = None,
                        cas_col: Optional[str] = None,
                        name_col: Optional[str] = None,
                        target_col: Optional[str] = None,
                        unit_col: Optional[str] = None,
                        inchikey_col: Optional[str] = None,
                        id_col: Optional[str] = None,
                        interactive: bool = True,
                        **kwargs) -> pd.DataFrame:
        """Load data from file"""
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File does not exist: {file_path}")

        loader_kwargs = kwargs.copy()
        loader_kwargs.update({
            'smiles_col': smiles_col,
            'cas_col': cas_col,
            'name_col': name_col,
            'target_col': target_col,
            'unit_col': unit_col,
            'inchikey_col': inchikey_col,
            'id_col': id_col
        })
        df = None
        suffix = os.path.splitext(file_path)[1].lower()
        if suffix == '.csv':
            df = pd.read_csv(file_path, **kwargs)
        elif suffix in {'.xls', '.xlsx'}:
            if 'sheet_name' in kwargs:
                df = pd.read_excel(file_path, **kwargs)
            else:
                with pd.ExcelFile(file_path) as xls:
                    sheet_names = xls.sheet_names
                if len(sheet_names) == 1:
                    df = pd.read_excel(file_path, sheet_name=sheet_names[0], **kwargs)
                else:
                    if not interactive:
                        raise ValueError(
                            "Multiple sheets found. Provide sheet_name. "
                            f"Available sheets: {sheet_names}"
                        )
                    print("The file contains multiple sheets. Please select one:")
                    for i, sheet_name in enumerate(sheet_names):
                        print(f"{i + 1}: {sheet_name}")

                    while True:
                        user_input = input("Enter the sheet number or name: ").strip()
                        try:
                            choice = int(user_input)
                            if 1 <= choice <= len(sheet_names):
                                df = pd.read_excel(file_path, sheet_name=sheet_names[choice - 1], **kwargs)
                                break
                            else:
                                print(f"Invalid number. Please enter a number between 1 and {len(sheet_names)}.")
                        except ValueError:
                            if user_input in sheet_names:
                                df = pd.read_excel(file_path, sheet_name=user_input)
                                break
                            else:
                                print(f"Invalid sheet name. Please enter a valid sheet name or number.")
        elif suffix == '.txt':
            df = pd.read_csv(file_path, sep='\t', **kwargs)
        elif suffix in {'.sdf', '.mol'}:
            df = DataHandler._load_sdf(file_path=file_path, **loader_kwargs)
        elif suffix == '.smi':
            df = DataHandler._load_smi(file_path=file_path, **loader_kwargs)
        else:
            logger.error("Only the .csv/.xlsx/.xls/.txt/.sdf/.mol/.smi file format is supported")
            raise ValueError("Unsupported file format")

        if df is None:
            raise ValueError("Failed to load data into DataFrame.")

        for col in df.select_dtypes(include=['object']).columns:
            numeric_col = pd.to_numeric(df[col], errors='coerce')
            if numeric_col.isnull().sum() > 0 and df[col].notnull().sum() > 0:
                df[col] = df[col].astype("string")

        check_map = {
            'SMILES': smiles_col, 'CAS': cas_col, 'Name': name_col, 'ID': id_col,
            'Target': target_col, 'Unit': unit_col, 'InChIKey': inchikey_col
        }
        for label, col_name in check_map.items():
            if col_name and col_name not in df.columns:
                if file_path.endswith('.sdf') and label == 'SMILES':
                    continue
                logger.error(f"{label} column '{col_name}' does not exist in the file. Available: {list(df.columns)}")
                raise KeyError(f"Missing column: {col_name}")

        return df

    @staticmethod
    def _load_from_list(smiles_list: List[str],
                        smiles_col: str) -> pd.DataFrame:
        """Load data from a list directly"""
        if not all(isinstance(s, str) for s in smiles_list):
            logger.error("The SMILES list must all be of string type")

        if smiles_col is None:
            logger.warning("The SMILES column name are not recommended to be left blank")
        data = {smiles_col: smiles_list}

        return pd.DataFrame(data)

    @staticmethod
    def _load_from_dataframe(df: pd.DataFrame,
                             smiles_col: str,
                             target_col: Optional[str] = None) -> pd.DataFrame:
        """Load data from an existing DataFrame"""
        required_cols = [smiles_col]
        if target_col:
            required_cols.append(target_col)

        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            logger.error(f"Missing necessary columns: {missing_cols}")

        return df.copy()

    @staticmethod
    def _load_sdf(file_path: str, smiles_col: Optional[str] = None, **kwargs) -> pd.DataFrame:
        from rdkit import Chem
        effective_smiles_col = smiles_col if smiles_col else 'smiles'

        def parse_mol_supplier(supplier):
            data_rows = []
            for record_number, mol in enumerate(supplier, 1):
                if mol is None:
                    data_rows.append({effective_smiles_col: None,
                                      'Input Record': record_number,
                                      'Structure Reliability': 'Invalid',
                                      'Import Error': 'SDF record could not be parsed or sanitized'})
                    continue
                try:
                    props = mol.GetPropsAsDict()
                    smi = Chem.MolToSmiles(mol, isomericSmiles=True)
                    declared = {
                        name: mol.GetProp(name) for name in mol.GetPropNames()
                        if name == effective_smiles_col
                        or name.strip().casefold() in DataHandler.SDF_SMILES_PROPERTIES
                    }
                    props['SDF Declared SMILES'] = json.dumps(declared, ensure_ascii=False, sort_keys=True)
                    props[effective_smiles_col] = smi
                    props['Input Record'] = record_number
                    problems = []
                    if mol.GetNumAtoms() == 0:
                        props['Structure Reliability'] = 'Invalid'
                        problems.append('SDF record contains no atoms')
                    else:
                        block_key = DataHandler._sdf_structure_key(mol)
                        parser = Chem.SmilesParserParams()
                        parser.removeHs = False
                        for name, value in declared.items():
                            if not value.strip():
                                continue
                            declared_mol = Chem.MolFromSmiles(value, parser)
                            if declared_mol is None or declared_mol.GetNumAtoms() == 0:
                                problems.append(f"SDF property {name!r} contains invalid SMILES {value!r}; molblock identity is {block_key!r}")
                            else:
                                declared_key = DataHandler._sdf_structure_key(declared_mol)
                                if declared_key != block_key:
                                    problems.append(f"SDF property {name!r} conflicts with molblock: declared {declared_key!r}; molblock {block_key!r}")
                        props['Structure Reliability'] = 'Unreliable' if problems else 'Reliable'
                    props['Import Error'] = '; '.join(problems) if problems else pd.NA
                    data_rows.append(props)
                except Exception as exc:
                    data_rows.append({effective_smiles_col: None,
                                      'Input Record': record_number,
                                      'Structure Reliability': 'Invalid',
                                      'Import Error': f'SDF record processing failed: {exc}'})
            return pd.DataFrame(data_rows)

        df = None
        try:
            with open(file_path, 'rb') as f:
                suppl = Chem.ForwardSDMolSupplier(f, removeHs=False, sanitize=True)
                df = parse_mol_supplier(suppl)
        except Exception as e:
            logger.error(f"Binary read failed: {str(e)}")
            try:
                suppl = Chem.SDMolSupplier(file_path, removeHs=False, sanitize=True)
                df = parse_mol_supplier(suppl)
            except Exception as final_e:
                logger.error(f"All parsing attempts failed: {str(final_e)}")
                raise

        if df is None or df.empty:
            raise ValueError(
                "Failed to parse molecules from the SDF file. The file might be corrupted or encoded strangely.")

        return df

    @staticmethod
    def _load_smi(file_path: str,
                  smiles_col: str = 'smiles',
                  id_col: Optional[str] = None,
                  **kwargs) -> pd.DataFrame:
        diptox_args = ['cas_col', 'name_col', 'target_col', 'unit_col', 'inchikey_col']
        for arg in diptox_args:
            kwargs.pop(arg, None)

        smiles_pos = kwargs.pop('smiles_pos', None)
        id_pos = kwargs.pop('id_pos', None)

        try:
            with open(file_path, encoding=kwargs.get('encoding', 'utf-8')) as f:
                first_line = f.readline()
            sep = kwargs.pop('sep', None) or ('\t' if '\t' in first_line else r'\s+')

            if 'header' in kwargs:
                header_infer = kwargs.pop('header')
            elif smiles_col and smiles_col in first_line:
                header_infer = 0
            else:
                header_infer = None

            df = pd.read_csv(file_path, sep=sep, header=header_infer, **kwargs)

            if header_infer is None:
                current_cols_str = [str(c) for c in df.columns]

                should_rename_smiles = True
                if smiles_col and str(smiles_col) in current_cols_str:
                    should_rename_smiles = False

                should_rename_id = True
                if id_col and str(id_col) in current_cols_str:
                    should_rename_id = False

                new_columns = list(df.columns)

                if smiles_pos is not None and isinstance(smiles_pos, int) and smiles_pos < len(new_columns):
                    smiles_idx = smiles_pos
                else:
                    smiles_idx = 0
                    if len(new_columns) >= 1:
                        from rdkit import Chem
                        from rdkit import RDLogger

                        RDLogger.DisableLog('rdApp.*')
                        for i in range(len(new_columns)):
                            val = str(df.iloc[0, i]).strip()
                            if not val or len(val) < 1:
                                continue
                            if Chem.MolFromSmiles(val) is not None:
                                smiles_idx = i
                                break
                        RDLogger.EnableLog('rdApp.*')

                if id_pos is not None and isinstance(id_pos, int) and id_pos < len(new_columns):
                    id_idx = id_pos
                else:
                    id_idx = 1 if smiles_idx == 0 else 0

                if should_rename_smiles and len(new_columns) > smiles_idx and smiles_col:
                    new_columns[smiles_idx] = smiles_col

                if should_rename_id and len(new_columns) > id_idx and id_col:
                    new_columns[id_idx] = id_col

                df.columns = new_columns
            df.columns = [str(c) for c in df.columns]

            if smiles_col and str(smiles_col) in df.columns:
                df[str(smiles_col)] = df[str(smiles_col)].astype(str).str.strip()
            return df

        except Exception as e:
            logger.error(f"Parsing the SMI file failed: {str(e)}")
            raise

    @staticmethod
    def save_data(df: pd.DataFrame, output_path: str, columns: list, smiles_col: str,
                  id_col: Optional[str] = None, *, interactive: bool = True,
                  source_smiles_col: Optional[str] = None) -> Optional[str]:
        """Save results, synchronizing exported source declarations with the SDF graph.

        ``source_smiles_col`` identifies a custom mapped source declaration when
        ``smiles_col`` selects a derived structure. Original values are archived.
        """
        output_path = os.fspath(output_path)
        suffix = os.path.splitext(output_path)[1].lower()
        if not interactive and suffix not in {'.csv', '.xls', '.xlsx', '.txt', '.sdf', '.smi'}:
            raise ValueError("Unsupported output file format. Use CSV, XLS, XLSX, TXT, SDF, or SMI.")
        missing_cols = [col for col in columns if col not in df.columns]
        if missing_cols:
            logger.error(f"Columns {missing_cols} not found in DataFrame")
            if not interactive:
                raise KeyError(f"Columns {missing_cols} not found in DataFrame")

        directory = os.path.dirname(output_path)
        if directory:
            os.makedirs(directory, exist_ok=True)

        while True:
            try:
                if suffix == '.csv':
                    from .unit_labels import excel_unit_label
                    export_frame = df[columns].apply(lambda column: column.map(excel_unit_label))
                    export_frame.to_csv(output_path, index=False, encoding='utf-8-sig')
                elif suffix in {'.xls', '.xlsx'}:
                    # Explicit text cells prevent Excel interpreting strings as formulas.
                    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
                        df[columns].to_excel(writer, index=False)
                        for row in writer.sheets['Sheet1'].iter_rows():
                            for cell in row:
                                if isinstance(cell.value, str):
                                    cell.data_type = 's'
                elif suffix == '.txt':
                    from .unit_labels import excel_unit_label
                    export_frame = df[columns].apply(lambda column: column.map(excel_unit_label))
                    export_frame.to_csv(output_path, index=False, sep='\t', encoding='utf-8-sig')
                elif suffix == '.sdf':
                    from rdkit import Chem
                    from rdkit.Chem import PandasTools
                    export_df = df.copy()
                    parser = Chem.SmilesParserParams()
                    parser.removeHs = False
                    molecules, statuses = [], []
                    declarations = {
                        name: [] for name in columns
                        if isinstance(name, str) and (
                            name == source_smiles_col
                            or name.strip().casefold() in DataHandler.SDF_SMILES_PROPERTIES
                        )
                    }
                    source_archives = []
                    archived_sources = False
                    for _, row in export_df.iterrows():
                        smiles = row[smiles_col]
                        mol = Chem.MolFromSmiles(smiles, parser) if isinstance(smiles, str) and smiles.strip() else None
                        if mol is None or mol.GetNumAtoms() == 0:
                            molecules.append(Chem.Mol())
                            statuses.append('Invalid structure placeholder')
                        else:
                            molecules.append(mol)
                            statuses.append('Written')
                        changed = {}
                        selected_key = DataHandler._sdf_structure_key(mol) if mol is not None and mol.GetNumAtoms() else None
                        for name, values in declarations.items():
                            value = row[name]
                            missing_value = pd.api.types.is_scalar(value) and pd.isna(value)
                            if selected_key is not None and not missing_value and str(value).strip():
                                declared_mol = Chem.MolFromSmiles(value, parser) if isinstance(value, str) else None
                                declared_key = DataHandler._sdf_structure_key(declared_mol) if declared_mol is not None else None
                                if declared_key != selected_key:
                                    changed[name] = str(value)
                                    value = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
                            values.append(value)
                        source_archive = row.get('SDF Source SMILES', pd.NA)
                        original_declarations = []
                        declared_json = row.get('SDF Declared SMILES')
                        if isinstance(declared_json, str):
                            try:
                                declared_sources = json.loads(declared_json)
                            except (TypeError, ValueError):
                                declared_sources = None
                            if isinstance(declared_sources, dict):
                                original_declarations = [
                                    (name, value) for name, value in declared_sources.items()
                                    if isinstance(name, str) and isinstance(value, str) and value.strip()
                                ]
                        source_values = [*original_declarations, *changed.items()]
                        if source_values:
                            archived_sources = True
                            archive = {}
                            missing_archive = pd.api.types.is_scalar(source_archive) and pd.isna(source_archive)
                            if not missing_archive and str(source_archive).strip():
                                try:
                                    archive = json.loads(str(source_archive))
                                except (TypeError, ValueError):
                                    archive = {'_previous': str(source_archive)}
                                if not isinstance(archive, dict):
                                    archive = {'_previous': archive}
                            for name, value in source_values:
                                if name not in archive:
                                    archive[name] = value
                                elif archive[name] != value:
                                    previous = archive[name] if isinstance(archive[name], list) else [archive[name]]
                                    archive[name] = previous if value in previous else [*previous, value]
                            source_archive = json.dumps(archive, ensure_ascii=False, sort_keys=True)
                        source_archives.append(source_archive)
                    export_df['ROMol'] = molecules
                    export_df['SDF Export Status'] = statuses
                    for name, values in declarations.items():
                        export_df[name] = values
                    if archived_sources:
                        export_df['SDF Source SMILES'] = source_archives
                    other = [c for c in columns if c != 'ROMol']
                    if 'SDF Export Status' not in other:
                        other.append('SDF Export Status')
                    if archived_sources and 'SDF Source SMILES' not in other:
                        other.append('SDF Source SMILES')
                    descriptor, temporary = tempfile.mkstemp(
                        prefix='.diptox-', suffix='.sdf', dir=os.path.abspath(directory or '.')
                    )
                    os.close(descriptor)
                    try:
                        PandasTools.WriteSDF(export_df[['ROMol'] + other], temporary,
                                             molColName='ROMol', properties=other)
                        os.replace(temporary, output_path)
                    finally:
                        if os.path.exists(temporary):
                            os.remove(temporary)
                elif suffix == '.smi':
                    if id_col is None or id_col not in columns or id_col not in df:
                        df[smiles_col].to_csv(output_path, sep='\t', header=True, index=False)
                    else:
                        df[[smiles_col, id_col]].to_csv(output_path, sep='\t', header=True, index=False)
                else:
                    logger.warning(f"Unsupported file format. The file will be saved as csv by default.")
                    output_path += '.csv'
                    suffix = '.csv'
                    df[columns].to_csv(output_path, index=False, encoding='utf-8')
                logger.info(f"File saved successfully: {output_path}")
                return output_path
            except (PermissionError, IOError, OSError) as e:
                logger.error(f"Unable to save file: {str(e)}")
                if not interactive:
                    raise
                choice = input("Do you want to save again? (Y/N): ").strip().lower()
                if choice in {'y', 'yes'}:
                    continue
                else:
                    logger.warning("The user canceled the save operation")
                    break
            except Exception as e:
                logger.error(f"Unknown error: {str(e)}")
                if not interactive:
                    raise
                break
