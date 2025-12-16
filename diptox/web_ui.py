# diptox/web_ui.py
import tempfile
from pathlib import Path
USER_TEMP_DIR = Path.home() / ".diptox" / "temp"
USER_TEMP_DIR.mkdir(parents=True, exist_ok=True)

import os
os.environ["DIPTOX_GUI_MODE"] = "true"

import streamlit as st
import pandas as pd
import threading
import time
from streamlit.runtime.scriptrunner import add_script_run_ctx, get_script_run_ctx
from diptox.core import DiptoxPipeline
from diptox import user_reg
from diptox.unit_processor import UnitProcessor

st.set_page_config(
    page_title="DiPTox GUI",
    page_icon="🧪",
    layout="wide"
)


def inject_custom_css():
    st.markdown(
        """
        <style>
        .stMultiSelect [data-baseweb="tag"] {
            background-color: #0078D4 !important; 
            color: white !important;
        }
        .stCheckbox input:checked + div[data-baseweb="checkbox"] {
            background-color: #0078D4 !important;
            border-color: #0078D4 !important;
        }
        .stSelectbox div[data-baseweb="select"] div[aria-selected="true"] {
            color: #0078D4;
        }
        .stAppDeployButton {
            visibility: hidden;
        }
        </style>
        """,
        unsafe_allow_html=True
    )


inject_custom_css()

# --- Session State ---
if 'pipeline' not in st.session_state:
    st.session_state.pipeline = DiptoxPipeline()
if 'df_preview' not in st.session_state:
    st.session_state.df_preview = None
if 'rule_logs' not in st.session_state:
    st.session_state.rule_logs = []
if 'user_cols' not in st.session_state:
    st.session_state.user_cols = set()
if 'gen_cols' not in st.session_state:
    st.session_state.gen_cols = set()
if 'export_data' not in st.session_state:
    st.session_state.export_data = None
if 'export_fname' not in st.session_state:
    st.session_state.export_fname = None
if 'dedup_condition_cols' not in st.session_state:
    st.session_state.dedup_condition_cols = []

pipeline = st.session_state.pipeline


def track_user_cols(cols):
    if not cols:
        return
    if isinstance(cols, str):
        cols = [cols]
    st.session_state.user_cols.update([c for c in cols if c])


def track_gen_cols(cols):
    if not cols:
        return
    if isinstance(cols, str):
        cols = [cols]
    st.session_state.gen_cols.update([c for c in cols if c])


# --- Side bar navigation ---
st.sidebar.title("DiPTox Control Panel")
step = st.sidebar.radio("Go to Step:",
                        ["Data Loading", "Preprocessing", "Web Requests", "Unit Standardization",
                         "Deduplication", "Search & Filter", "Export"])


def update_preview():
    if pipeline.df is not None:
        st.session_state.df_preview = pipeline.df
    else:
        st.session_state.df_preview = None


def add_rule_log(message):
    st.session_state.rule_logs.append(message)
    st.toast(message, icon="✅")


if not user_reg.is_registered_or_skipped():
    with st.expander("👋 DiPTox Community Check-in (Optional)", expanded=True):
        st.markdown("""
        <small>Hi! To help us understand our user base and improve DiPTox, we would appreciate it if you could share who you are.
        This helps us with academic impact assessment. This is strictly optional.</small>
        """, unsafe_allow_html=True)

        with st.form("user_reg_form"):
            col_u1, col_u2, col_u3 = st.columns(3)
            reg_name = col_u1.text_input("Name")
            reg_aff = col_u2.text_input("Affiliation")
            reg_email = col_u3.text_input("Email")

            # Button row
            col_btn1, col_btn2 = st.columns([1, 5])
            is_submit = col_btn1.form_submit_button("🚀 Submit", type="primary")

        # Skip button outside the form
        if st.button("Skip (Don't ask again)"):
            user_reg.save_status("skipped")
            st.rerun()

        if is_submit:
            if not reg_name or not reg_aff:
                st.error("Please fill in at least Name and Affiliation.")
            else:
                with st.spinner("Sending..."):
                    success, msg = user_reg.submit_info(reg_name, reg_aff, reg_email)
                    if success:
                        st.success(msg)
                        import time

                        time.sleep(1.5)
                        st.rerun()
                    else:
                        st.warning(f"Note: {msg}")
                        # Even if it fails, mark as registered to avoid annoying the user
                        user_reg.save_status("registered_offline")
                        time.sleep(2)
                        st.rerun()
    st.divider()


# Data Loading
if step == "Data Loading":
    st.header("Data Loading")

    tab_file, tab_text = st.tabs(["Upload File", "Paste SMILES List"])

    with tab_file:
        st.markdown("Supports: `.csv`, `.xlsx`, `.xls`, `.txt`, `.sdf`, `.smi`, `.mol`")
        st.warning("⚠️ Note: For files larger than 200MB, please use the Python script directly.")

        uploaded_file = st.file_uploader("Upload Dataset", type=['csv', 'xlsx', 'xls', 'txt', 'smi', 'sdf', 'mol'])

        if uploaded_file:
            current_pid = os.getpid()
            file_id = f"{current_pid}_{uploaded_file.name}_{uploaded_file.size}"
            suffix = os.path.splitext(uploaded_file.name)[1].lower()

            if "current_file_id" not in st.session_state:
                st.session_state["current_file_id"] = None
            if "current_temp_path" not in st.session_state:
                st.session_state["current_temp_path"] = None

            if st.session_state["current_file_id"] != file_id:
                old_path = st.session_state["current_temp_path"]
                if old_path and os.path.exists(old_path):
                    try:
                        os.remove(old_path)
                    except:
                        pass
                with tempfile.NamedTemporaryFile(
                    dir=str(USER_TEMP_DIR),
                    delete=False,
                    suffix=suffix,
                    prefix=f"tmp_{current_pid}_"
                ) as tmp_file:
                    tmp_file.write(uploaded_file.getbuffer())
                    st.session_state["current_temp_path"] = tmp_file.name
                st.session_state["current_file_id"] = file_id

            temp_filename = st.session_state["current_temp_path"]
            if temp_filename and os.path.exists(temp_filename):
                file_ext = suffix
                detected_cols = []
                selected_sheet = None
                has_header = True
                with st.expander("Column Configuration", expanded=True):
                    if file_ext in ['.sdf', '.mol']:
                        st.info(
                            f"💡 **{file_ext.upper()} File Detected**: You can leave the **'SMILES Column'** as **(None)**. "
                            "Molecules will be automatically parsed from the structure block.")
                    if file_ext in ['.csv', '.txt', '.smi']:
                        has_header = st.checkbox("File has a header row?", value=True,
                                                 help="Uncheck this if the first row contains data, not column names.")
                    try:
                        # === 1. Excel ===
                        if file_ext in ['.xlsx', '.xls']:
                            with pd.ExcelFile(temp_filename) as xl:
                                sheet_names = xl.sheet_names
                                if len(sheet_names) > 1:
                                    st.info(f"Multi-sheet Excel detected.")
                                    selected_sheet = st.selectbox("Select Excel Sheet:", sheet_names)
                                else:
                                    selected_sheet = sheet_names[0]
                                df_preview = pd.read_excel(temp_filename, sheet_name=selected_sheet, nrows=5)
                                detected_cols = list(df_preview.columns)
                        # === 2. SDF ===
                        elif file_ext in ['.sdf', '.mol']:
                            from rdkit import Chem
                            suppl = None
                            try:
                                suppl = Chem.ForwardSDMolSupplier(open(temp_filename, 'rb'))
                                mol = next(suppl) if suppl else None
                            except:
                                mol = None
                            if mol:
                                detected_cols = list(mol.GetPropsAsDict().keys())
                            else:
                                st.warning("Could not preview SDF properties (File might be empty or encrypted).")
                                detected_cols = []
                        # === 3. CSV/TXT/SMI ===
                        else:
                            header_arg = 0 if has_header else None
                            sep = '\t' if file_ext in ['.txt', '.smi'] else ','
                            try:
                                df_preview = pd.read_csv(temp_filename, sep=sep, nrows=5, header=header_arg)
                                detected_cols = list(df_preview.columns)
                            except:
                                try:
                                    df_preview = pd.read_csv(temp_filename, nrows=5, header=header_arg)
                                    detected_cols = list(df_preview.columns)
                                except:
                                    detected_cols = []
                        detected_cols = [str(c) for c in detected_cols]
                    except Exception as e:
                        st.error(f"Failed to parse headers: {e}")
                        detected_cols = []
                    if detected_cols:
                        st.success(f"✅ Detected {len(detected_cols)} columns.")
                        options_with_none = ["(None)"] + detected_cols
                        c1, c2 = st.columns(2)
                        # 1. SMILES Column
                        smiles_col = c1.selectbox("SMILES Column (Optional)", options_with_none, index=0)
                        # 2. Target Column
                        name_col = c2.selectbox("Name Column (Optional)", options_with_none, index=0)
                        c3, c4 = st.columns(2)
                        # 3. CAS Column
                        cas_col = c3.selectbox("CAS Column (Optional)", options_with_none, index=0)
                        # 4. ID Column
                        id_col = c4.selectbox("ID Column (Optional)", options_with_none, index=0)
                        c5, c6 = st.columns(2)
                        # 5. Name Column
                        target_col = c5.selectbox("Target Column (Optional)", options_with_none, index=0)
                        # 6. Unit Column
                        unit_col = c6.selectbox("Unit Column (Optional)", options_with_none, index=0)
                    else:
                        st.warning("Could not detect columns automatically. Please enter manually.")
                        with st.expander("Column Configuration", expanded=True):
                            c1, c2 = st.columns(2)
                            smiles_col = c1.text_input("SMILES Column (Optional)", value="")
                            name_col = c2.text_input("Name Column (Optional)", value="")
                            c3, c4 = st.columns(2)
                            cas_col = c3.text_input("CAS Column (Optional)", value="")
                            id_col = c4.text_input("ID Column (for .smi)", value="")
                            c5, c6 = st.columns(2)
                            target_col = c5.text_input("Target Column (Optional)", value="")
                            unit_col = c6.text_input("Unit Column (Optional)", value="")
                    current_selections = {
                        "SMILES Column": smiles_col, "Name Column": name_col, "CAS Column": cas_col,
                        "ID Column": id_col, "Target Column": target_col, "Unit Column": unit_col
                    }
                    active_selections = {}
                    for label, val in current_selections.items():
                        if val and val != "(None)" and str(val).strip() != "":
                            active_selections[label] = val
                    seen_cols = {}
                    dupe_warnings = []
                    for label, col_name in active_selections.items():
                        if col_name in seen_cols:
                            prev_label = seen_cols[col_name]
                            dupe_warnings.append(
                                f"• **{prev_label}** and **{label}** are both mapped to column: `{col_name}`")
                        else:
                            seen_cols[col_name] = label
                    if dupe_warnings:
                        st.warning("⚠️ **Duplicate Column Assignment Detected**:\n" + "\n".join(dupe_warnings))
            if st.button("Load Data from File", type="primary"):
                try:
                    kwargs = {}
                    if selected_sheet:
                        kwargs['sheet_name'] = selected_sheet
                    if file_ext in ['.csv', '.txt', '.smi']:
                        kwargs['header'] = 0 if has_header else None
                    def clean_col_name(val):
                        if val == "(None)":
                            return None
                        if isinstance(val, str) and val.strip() == "":
                            return None
                        return val
                    final_smiles = clean_col_name(smiles_col)
                    final_target = clean_col_name(target_col)
                    final_unit = clean_col_name(unit_col)
                    final_cas = clean_col_name(cas_col)
                    final_name = clean_col_name(name_col)
                    final_id = clean_col_name(id_col)
                    pipeline.load_data(
                        input_data=temp_filename,
                        smiles_col=final_smiles,
                        cas_col=final_cas,
                        target_col=final_target,
                        unit_col=final_unit,
                        name_col=final_name,
                        id_col=final_id,
                        **kwargs
                    )
                    cols_to_track = [
                        pipeline.smiles_col,
                        pipeline.target_col,
                        pipeline.unit_col,
                        pipeline.cas_col,
                        pipeline.name_col,
                        pipeline.id_col
                    ]
                    track_user_cols(cols_to_track)
                    update_preview()
                    st.success(f"Successfully loaded {len(pipeline.df)} records!")
                except Exception as e:
                    st.error(f"Error loading data: {str(e)}")

    with tab_text:
        st.info("Paste a list of SMILES strings, one per line.")

        raw_text = st.text_area("SMILES List", height=300, placeholder="C\nCC\nCCC\nCN(C)C=O")

        c_t1, c_t2 = st.columns([1, 1])
        custom_smiles_col = c_t1.text_input("Resulting Column Name", value="Smiles")

        if st.button("Load Data from Text", type="primary"):
            if not raw_text.strip():
                st.error("Please enter some SMILES strings.")
            else:
                try:
                    smiles_list = [line.strip() for line in raw_text.split('\n') if line.strip()]

                    if not smiles_list:
                        st.error("No valid lines found.")
                    else:
                        pipeline.load_data(
                            input_data=smiles_list,
                            smiles_col=custom_smiles_col
                        )

                        track_user_cols(custom_smiles_col)

                        update_preview()
                        st.success(f"Successfully loaded {len(pipeline.df)} SMILES!")

                except Exception as e:
                    st.error(f"Error loading text data: {str(e)}")

# Preprocessing (Rules & Run)
elif step == "Preprocessing":
    st.header("Preprocessing & Standardization")

    if pipeline.df is None:
        st.warning("Please load data first.")
    else:
        st.markdown("### A. Rule Management (Optional)")

        with st.expander("🛠️ Manage & View Chemical Rules", expanded=False):
            st.caption("Add or remove rules. Changes are applied immediately.")

            # --- 1. Atoms ---
            c1, c2, c3 = st.columns([3, 1, 1], vertical_alignment="bottom")
            atom_input = c1.text_input("Valid Atoms", placeholder="e.g., Si, Zr", key="atom_in")

            if c2.button("➕ Add", key="btn_add_atom", width="stretch"):
                if not atom_input.strip():
                    st.error("Empty input")
                else:
                    atoms = [x.strip() for x in atom_input.split(',') if x.strip()]
                    failed = pipeline.manage_atom_rules(atoms=atoms, add=True)
                    if failed:
                        st.error(f"Failed to add invalid atoms: {failed}")
                    else:
                        add_rule_log(f"Added atom(s): {atom_input}")
                        st.rerun()

            if c3.button("➖ Del", key="btn_del_atom", width="stretch"):
                if not atom_input.strip():
                    st.error("Empty input")
                else:
                    atoms = [x.strip() for x in atom_input.split(',') if x.strip()]
                    failed = pipeline.manage_atom_rules(atoms=atoms, add=False)
                    if failed:
                        st.warning(f"Could not find/remove: {failed}")
                    else:
                        add_rule_log(f"Removed atom(s): {atom_input}")
                        st.rerun()

            # --- 2. Salts ---
            c1, c2, c3 = st.columns([3, 1, 1], vertical_alignment="bottom")
            salt_input = c1.text_input("Custom Salts (SMARTS)", placeholder="e.g., [Hg+2]", key="salt_in")

            if c2.button("➕ Add", key="btn_add_salt", width="stretch"):
                if not salt_input.strip():
                    st.error("Empty input")
                else:
                    salts = [x.strip() for x in salt_input.split(',') if x.strip()]
                    failed = pipeline.manage_default_salt(salts=salts, add=True)
                    if failed:
                        st.error(f"Invalid SMARTS patterns: {failed}")
                    else:
                        add_rule_log(f"Added salt: {salt_input}")
                        st.rerun()

            if c3.button("➖ Del", key="btn_del_salt", width="stretch"):
                if not salt_input.strip():
                    st.error("Empty input")
                else:
                    salts = [x.strip() for x in salt_input.split(',') if x.strip()]
                    pipeline.manage_default_salt(salts=salts, add=False)
                    add_rule_log(f"Removed salt: {salt_input}")
                    st.rerun()

            # --- 3. Solvents ---
            c1, c2, c3 = st.columns([3, 1, 1], vertical_alignment="bottom")
            solv_input = c1.text_input("Custom Solvents (SMILES)", placeholder="e.g., CCC", key="solv_in")

            if c2.button("➕ Add", key="btn_add_solv", width="stretch"):
                if not solv_input.strip():
                    st.error("Empty input")
                else:
                    solvs = [x.strip() for x in solv_input.split(',') if x.strip()]
                    failed = pipeline.manage_default_solvent(solvents=solvs, add=True)
                    if failed:
                        st.error(f"Invalid SMILES: {failed}")
                    else:
                        add_rule_log(f"Added solvent: {solv_input}")
                        st.rerun()

            if c3.button("➖ Del", key="btn_del_solv", width="stretch"):
                if not solv_input.strip():
                    st.error("Empty input")
                else:
                    solvs = [x.strip() for x in solv_input.split(',') if x.strip()]
                    pipeline.manage_default_solvent(solvents=solvs, add=False)
                    add_rule_log(f"Removed solvent: {solv_input}")
                    st.rerun()

            # --- 4. Neutralization ---
            st.markdown("**Neutralization Rules**")
            c1, c2, c3, c4 = st.columns([1.5, 1.5, 0.5, 0.5], vertical_alignment="bottom")
            react = c1.text_input("Reactant (SMARTS)", key="neu_r")
            prod = c2.text_input("Product (SMILES)", key="neu_p")

            if c3.button("➕", key="btn_add_neu", help="Add Rule", width="stretch"):
                success = pipeline.add_neutralization_rule(react, prod)
                if success:
                    add_rule_log(f"Added rule: {react} -> {prod}")
                    st.rerun()
                else:
                    st.error("Invalid Rule (Check SMARTS/SMILES)")

            if c4.button("➖", key="btn_del_neu", help="Remove Rule", width="stretch"):
                success = pipeline.remove_neutralization_rule(react)
                if success:
                    add_rule_log(f"Removed rule: {react}")
                    st.rerun()
                else:
                    st.warning(f"Rule not found: {react}")

            # --- Rule Viewer Tabs ---
            st.markdown("---")
            current_rules = pipeline.chem_processor.get_current_rules_dict()
            t1, t2, t3, t4 = st.tabs(["Atoms", "Salts", "Solvents", "Neutralization"])
            with t1:
                st.code(", ".join(current_rules['atoms']), language="text")
            with t2:
                st.dataframe(pd.DataFrame(current_rules['salts'], columns=["Salt (SMARTS)"]), hide_index=True,
                             width="stretch")
            with t3:
                st.dataframe(pd.DataFrame(current_rules['solvents'], columns=["Solvent (SMILES)"]), hide_index=True,
                             width="stretch")
            with t4:
                st.dataframe(pd.DataFrame(current_rules['neutralization'], columns=["Reactant", "Product"]),
                             hide_index=True, width="stretch")

        # --- B. Execute Preprocessing ---
        st.markdown("### B. Execute Preprocessing")

        col_clean_card, col_norm_card = st.columns(2)

        with col_clean_card:
            with st.container(border=True):
                st.markdown("#### Cleaning")
                st.caption("Removal of unwanted fragments")

                # Row 1: Salts | Solvents
                c_cl1, c_cl2 = st.columns(2)
                rm_salts = c_cl1.checkbox("Remove Salts", value=True, help="Strip salt fragments.")
                rm_solvs = c_cl2.checkbox("Remove Solvents", value=True, help="Strip solvent fragments.")

                # Row 2: Inorganic | Mixtures
                c_cl3, c_cl4 = st.columns(2)
                rm_inorg = c_cl3.checkbox("Remove Inorganic", value=True, help="Remove inorganic molecules.")
                rm_mixt = c_cl4.checkbox("Remove Mixtures", value=False,
                                         help="Handle molecules with multiple disconnected fragments.")

                # Row 3: Keep Largest | HAC Threshold
                c_cl5, c_cl6 = st.columns(2, vertical_alignment="center")

                # 3.1 Keep Largest
                keep_large = c_cl5.checkbox("Keep Largest", value=True, disabled=not rm_mixt,
                                            help="Keep the largest fragment.")

                # 3.2 HAC Config
                with c_cl6:
                    c_hac_chk, c_hac_inp = st.columns([1.5, 1], vertical_alignment="center")
                    use_hac = c_hac_chk.checkbox("HAC Filter", value=True, disabled=not rm_mixt,
                                                 help="Filter by Heavy Atom Count.")
                    hac_val = c_hac_inp.number_input("HAC", value=3, min_value=0, label_visibility="collapsed",
                                                     disabled=(not rm_mixt or not use_hac))

                final_hac_threshold = hac_val if (rm_mixt and use_hac) else 0

        with col_norm_card:
            with st.container(border=True):
                st.markdown("#### Normalization")
                st.caption("Standardization & Validation")

                # RDKit Sanitize
                c_n0, c_n = st.columns(2)
                sanitize = c_n0.checkbox("RDKit Sanitize", value=True,
                                         help="Ensure chemical validity using RDKit (valence, aromaticity, etc.).")

                # Row 2: Hs | Isotopes
                c_n1, c_n2 = st.columns(2)
                rm_hs = c_n1.checkbox("Remove Hs", value=True, help="Remove explicit hydrogen atoms.")
                rm_iso = c_n2.checkbox("Remove Isotopes", value=True,
                                       help="Remove isotopic information (e.g., [13C] -> C).")

                # Row 3: Stereo | Radicals
                c_n3, c_n4 = st.columns(2)
                rm_stereo = c_n3.checkbox("Remove Stereo", value=False,
                                          help="Remove stereochemistry information (e.g., @, /\\).")
                rej_rad = c_n4.checkbox("Reject Radicals", value=True,
                                        help="Discard molecules containing radical electrons.")

                # Row 4: Neutralize | Reject Non-Neutral
                c_neu1, c_neu2 = st.columns(2)
                neu = c_neu1.checkbox("Neutralize Charges", value=True,
                                      help="Neutralize charges based on predefined rules.")
                rej_non_neu = c_neu2.checkbox("Reject Non-Neutral", value=False, disabled=not neu,
                                              help="Discard if charge remains after adjustment.")

                # Row 5: Check Valid | Strict Mode
                c_chk1, c_chk2 = st.columns(2)
                chk_valid = c_chk1.checkbox("Check Valid Atoms", value=False,
                                            help="Validate against the allowed atom list.")
                strict = c_chk2.checkbox("Strict Mode", value=False, disabled=not chk_valid,
                                         help="Reject entire molecule if invalid atoms found.")

        # --- C. Post-Processing ---
        st.markdown("### C. Post-Processing")
        with st.container(border=True):
            calc_inchi = st.checkbox("Calculate InChI", value=False, help="Generate InChI strings after processing.")

        if st.button("🚀 Run Preprocessing", type="primary"):
            progress_container = st.empty()

            def update_progress(current, total):
                if total > 0:
                    with progress_container.container():
                        st.markdown("#### ⚙️ Processing molecules...")
                        st.progress(min(current / total, 1.0))
                        st.caption(f"Progress: {current}/{total}")

            try:
                update_progress(0, 100)
                pipeline.preprocess(
                    remove_salts=rm_salts,
                    remove_solvents=rm_solvs,
                    remove_mixtures=rm_mixt,
                    remove_inorganic=rm_inorg,
                    neutralize=neu,
                    reject_non_neutral=rej_non_neu,
                    check_valid_atoms=chk_valid,
                    strict_atom_check=strict,
                    remove_stereo=rm_stereo,
                    remove_isotopes=rm_iso,
                    remove_hs=rm_hs,
                    keep_largest_fragment=keep_large,
                    hac_threshold=final_hac_threshold,
                    sanitize=sanitize,
                    reject_radical_species=rej_rad,
                    progress_callback=update_progress
                )
                new_cols = ['Canonical SMILES', 'Is Valid', 'Processing Log']
                if calc_inchi:
                    pipeline.calculate_inchi()
                    new_cols.append('InChI')
                track_gen_cols(new_cols)
                update_preview()
                st.success("Preprocessing complete!")
            except Exception as e:
                st.error(f"Preprocessing failed: {str(e)}")
            finally:
                progress_container.empty()

# Web Requests
elif step == "Web Requests":
    st.header("Web Requests")

    if pipeline.df is None:
        st.warning("Please load data first.")
    else:
        st.markdown("Retrieve data from multiple sources simultaneously.")

        st.info(
            "💡 Tip: To **STOP** a running request, click the **'Stop'** button in the top-right corner of this webpage.")

        col_req1, col_req2 = st.columns(2)

        with col_req1:
            sources = st.multiselect(
                "Select Data Sources (Priority Order)",
                ['pubchem', 'chemspider', 'comptox', 'cas', 'cactus', 'chembl'],
                default=['pubchem']
            )
            if sources:
                st.caption(f"🔍 Source Priority: {' → '.join(sources)}")
            else:
                st.caption("⚠️ No source selected")

            req_types = st.multiselect(
                "Request Properties (Output)",
                ["smiles", "cas", "name", "iupac", "mw"],
                default=["cas", "iupac"]
            )

            # API Keys
            chemspider_key = None
            comptox_key = None
            cas_key = None

            if 'chemspider' in sources:
                chemspider_key = st.text_input("ChemSpider API Key", type="password")
            if 'comptox' in sources:
                comptox_key = st.text_input("CompTox API Key", type="password")
            if 'cas' in sources:
                cas_key = st.text_input("CAS Common Chemistry API Key", type="password")

        with col_req2:
            send_type = st.multiselect(
                "Send Identifier Input (Priority Order)",
                ["smiles", "cas", "name"],
                default=["smiles"]
            )
            if send_type:
                st.caption(f"📥 Input Priority: {' → '.join(send_type)}")
            else:
                st.caption("⚠️ No input type selected")

            max_workers = st.number_input("Max Workers (Threads)", min_value=1, max_value=16, value=4, step=1,
                                          help="Number of concurrent requests.")

        with st.expander("⚙️ Advanced Configuration (Timeouts, Retries, Limits)", expanded=False):
            st.caption("Fine-tune network behavior to handle API rate limits or unstable connections.")

            c_adv1, c_adv2, c_adv3 = st.columns(3)
            interval = c_adv1.number_input("Request Interval (s)", value=0.3, min_value=0.0, step=0.1,
                                           help="Time to wait between requests to avoid rate limits.")
            retries = c_adv2.number_input("Max Retries", value=3, min_value=0,
                                          help="Number of attempts if a request fails.")
            retry_delay = c_adv3.number_input("Retry Delay (s)", value=30, min_value=0,
                                              help="Time to wait before retrying a failed request.")

            c_adv4, c_adv5, c_adv6 = st.columns(3)
            batch_limit = c_adv4.number_input("Batch Limit", value=1500, min_value=0,
                                              help="Pause execution after X requests to rest (0 = no limit).")
            rest_duration = c_adv5.number_input("Rest Duration (s)", value=300, min_value=0,
                                                help="Time to pause (in seconds) when batch limit is reached.")
            mode_selection = c_adv6.selectbox(
                "API Mode Strategy",
                options=["Auto (SDK)", "Force API"],
                index=0,
                help="Auto: Use official SDKs (Faster/Stable).\nForce API: Bypass SDKs, use raw HTTP requests."
            )

            force_api = True if mode_selection == "Force API" else False

        if st.button("Start Web Request", type="primary"):
            if not sources:
                st.error("Please select at least one data source.")
            elif not send_type:
                st.error("Please select at least one input identifier type.")
            else:
                main_ctx = get_script_run_ctx()
                progress_container = st.empty()

                def update_progress(current, total):
                    if main_ctx:
                        add_script_run_ctx(threading.current_thread(), main_ctx)

                    if total > 0:
                        percent = min(current / total, 1.0)
                        with progress_container.container():
                            st.progress(percent)
                            st.caption(f"Querying batch: {current}/{total}")

                def update_status(message):
                    if main_ctx:
                        add_script_run_ctx(threading.current_thread(), main_ctx)

                    with progress_container.container():
                        st.warning(message, icon="⏳")

                try:
                    with st.spinner(f"Querying... (Check terminal for logs)"):
                        pipeline.config_web_request(
                            sources=sources,
                            chemspider_api_key=chemspider_key,
                            comptox_api_key=comptox_key,
                            cas_api_key=cas_key,
                            max_workers=max_workers,
                            interval=interval,
                            retries=retries,
                            delay=retry_delay,
                            batch_limit=batch_limit,
                            rest_duration=rest_duration,
                            force_api_mode=force_api,
                            status_callback=update_status
                        )
                        pipeline.web_request(send=send_type, request=req_types, progress_callback=update_progress)
                        new_web_cols = [f"{prop}_from_web" for prop in req_types]
                        new_web_cols.extend(['Query_Status', 'Data_Source', 'Query_Method'])
                        track_gen_cols(new_web_cols)
                        update_preview()
                        st.success("Web request finished.")
                except Exception as e:
                    st.error(f"Web request failed: {str(e)}")
                finally:
                    progress_container.empty()

# Unit Standardization
elif step == "Unit Standardization":
    st.header("Unit Standardization")

    if pipeline.df is None:
        st.warning("Please load data first.")
    elif not pipeline.target_col or not pipeline.unit_col:
        st.warning("This step requires both a 'Target Column' and a 'Unit Column' to be specified during Data Loading.")
    else:
        st.info(f"Processing Target Column: **{pipeline.target_col}** with units from: **{pipeline.unit_col}**")

        unique_units = [u for u in pipeline.df[pipeline.unit_col].dropna().unique() if u]

        if len(unique_units) <= 1:
            st.success("Only one unit was detected, or no units were found. This step is not required.")
        else:
            st.markdown("### Define Conversion Rules")
            st.caption(
                "Select a standard unit and provide formulas to convert others to it. Leave a formula blank to remove data with that unit.")

            st.info(f"Detected Units: **{', '.join(unique_units)}**")

            standard_unit = st.selectbox("Select the Standard Unit", options=unique_units)

            conversion_rules = {}
            temp_processor = UnitProcessor()

            with st.container(border=True):
                st.markdown("##### Conversion Formulas")

                other_units = [u for u in unique_units if u != standard_unit]
                col1, col2, col3 = st.columns(3)

                for i, unit in enumerate(other_units):
                    target_col = col1 if i % 3 == 0 else (col2 if i % 3 == 1 else col3)
                    default_rule = temp_processor.get_rule(unit, standard_unit)
                    pre_fill_value = default_rule if default_rule else ""
                    help_msg = "Built-in rule detected." if default_rule else "No built-in rule. Please enter manually."
                    with target_col:
                        formula = st.text_input(
                            label=f"**{unit}** → **{standard_unit}**",
                            value=pre_fill_value,
                            placeholder="e.g., x * 1000",
                            help=f"{help_msg} Use 'x' as value. Clear to skip.",
                            key=f"rule_{unit}_to_{standard_unit}"
                        )
                        if formula:
                            conversion_rules[(unit, standard_unit)] = formula

            if st.button("🚀 Run Unit Standardization", type="primary"):
                try:
                    with st.spinner("Applying unit conversions..."):
                        pipeline.standardize_units(
                            standard_unit=standard_unit,
                            conversion_rules=conversion_rules
                        )

                        new_target_col = pipeline.target_col

                        track_gen_cols([new_target_col])
                        update_preview()
                        st.success(f"Unit standardization complete! New column '{new_target_col}' created.")

                except ValueError as e:
                    st.error(f"Error: {e}. Please ensure all required conversion formulas are provided and are valid.")
                except Exception as e:
                    st.error(f"An unexpected error occurred: {e}")

# Deduplication
elif step == "Deduplication":
    st.header("Dataset Deduplication")

    if pipeline.df is None:
        st.warning("Please load data first.")
    else:
        avail_cols = list(pipeline.df.columns)
        c_cond_main, c_cond_opt = st.columns([3, 1])

        with c_cond_main:
            condition_cols = st.multiselect("Condition Columns (e.g., pH, Temperature)", avail_cols)

        with c_cond_opt:
            dropna_selection = st.selectbox(
                "Drop NaN Conditions",
                options=["False", "True"],
                index=0,
                help="Drop rows with missing condition values."
            )
            dropna_conditions = False if dropna_selection == "False" else True

        data_type = st.selectbox("Data Type", ["continuous", "discrete", "smiles"])

        method = "auto"
        priority_list = None
        log_transform = False

        if data_type == "continuous":
            c1, c2 = st.columns(2)
            with c1:
                method = st.selectbox("Deduplication Method", ["auto", "3sigma", "IQR"])
            with c2:
                mode_selection = st.selectbox(
                    "Apply -log10 Transformation",
                    options=["False", "True"],
                    index=0,
                    help="Transforms target values using `y = -log10(x)`. Rows with non-positive values will be removed."
                )
                log_transform = False if mode_selection == "False" else True
        elif data_type == "discrete":
            method = "vote"
            st.info("Priority Rule: If specified values exist in a duplicate group, they are selected first.")
            priority_input = st.text_area("Priority Values (comma separated, e.g., Active, Intermediate)")
            if priority_input:
                priority_list = [x.strip() for x in priority_input.split(",") if x.strip()]

        if st.button("Run Deduplication", type="primary"):
            progress_container = st.empty()

            def update_progress(current, total):
                if total > 0:
                    with progress_container.container():
                        st.progress(min(current / total, 1.0))
                        st.caption(f"Processing Groups: {current}/{total}")

            try:
                st.session_state.dedup_condition_cols = condition_cols if condition_cols else []
                pipeline.config_deduplicator(
                    condition_cols=condition_cols if condition_cols else None,
                    data_type=data_type,
                    method=method,
                    priority=priority_list,
                    log_transform=log_transform,
                    dropna_conditions=dropna_conditions
                )
                pipeline.dataset_deduplicate()
                final_target_col = pipeline.target_col
                cols_to_light = ['Deduplication Strategy']
                if final_target_col and final_target_col in pipeline.df.columns:
                    cols_to_light.append(final_target_col)
                track_gen_cols(cols_to_light)
                update_preview()
                st.success("Deduplication complete.")
            except ValueError as e:
                if "Unit standardization is required" in str(e):
                    st.error(
                        "Deduplication failed: Unit standardization is required but has not been completed. Please go to the 'Unit Standardization' step and run it first.")
                else:
                    st.error(f"Deduplication failed with an error: {e}")
            except Exception as e:
                st.error(f"Deduplication failed: {str(e)}")
            finally:
                progress_container.empty()

# Search & Filter
elif step == "Search & Filter":
    st.header("Search & Filter")

    if pipeline.df is None:
        st.warning("Please load data first.")
    else:
        tab_search, tab_filter = st.tabs(["Substructure Search", "Atom Count Filter"])

        with tab_search:
            st.markdown("Annotate molecules containing specific substructures.")
            query = st.text_input("Enter Pattern (SMARTS/SMILES)")
            is_smarts = st.checkbox("Is SMARTS Pattern?", value=True)

            if st.button("Search", type="primary"):
                try:
                    pipeline.substructure_search(query, is_smarts=is_smarts)
                    track_gen_cols(f'Substructure_{query}')
                    update_preview()
                    st.success(f"Search complete. Check for column 'Substructure_{query}'")
                except Exception as e:
                    st.error(f"Search failed: {str(e)}")

        with tab_filter:
            st.markdown("Filter dataset based on atom counts.")

            c_f1, c_f2 = st.columns(2)
            min_heavy = c_f1.number_input("Min Heavy Atoms", value=0)
            max_heavy = c_f2.number_input("Max Heavy Atoms", value=999)

            c_f3, c_f4 = st.columns(2)
            min_total = c_f3.number_input("Min Total Atoms", value=0)
            max_total = c_f4.number_input("Max Total Atoms", value=999)

            use_heavy = st.checkbox("Apply Heavy Atom Filter", value=True)
            use_total = st.checkbox("Apply Total Atom Filter", value=False)

            if st.button("Apply Filter", type="primary"):
                try:
                    kwargs = {}
                    if use_heavy:
                        kwargs['min_heavy_atoms'] = min_heavy
                        kwargs['max_heavy_atoms'] = max_heavy
                    if use_total:
                        kwargs['min_total_atoms'] = min_total
                        kwargs['max_total_atoms'] = max_total

                    if kwargs:
                        pipeline.filter_by_atom_count(**kwargs)
                        update_preview()
                        st.success(f"Filter applied. Remaining rows: {len(pipeline.df)}")
                    else:
                        st.warning("No filter criteria selected.")
                except Exception as e:
                    st.error(f"Filtering failed: {str(e)}")

# Export
elif step == "Export":
    st.header("Export Results")

    if pipeline.df is None:
        st.warning("No data to export.")
    else:
        st.markdown(f"**Current Dataset Shape:** {pipeline.df.shape[0]} rows, {pipeline.df.shape[1]} columns")

        col_file1, col_file2 = st.columns([1, 3])
        with col_file1:
            save_fmt = st.selectbox("File Format", ["csv", "xlsx", "txt", "sdf", "smi"])
        with col_file2:
            base_filename = st.text_input("Output Filename (Without extension)", value="diptox_processed")
        clean_base_name = os.path.splitext(base_filename)[0]
        filename = f"{clean_base_name}.{save_fmt}"
        st.info(f"💾 File will be saved as: **{filename}**")

        st.divider()
        st.subheader("Select Columns")
        all_cols = list(pipeline.df.columns)

        if 'export_selected_cols' not in st.session_state:
            st.session_state.export_selected_cols = all_cols

        def get_recommended_cols():
            rec = []
            if pipeline.smiles_col and pipeline.smiles_col in all_cols: rec.append(pipeline.smiles_col)
            if pipeline.id_col and pipeline.id_col in all_cols: rec.append(pipeline.id_col)
            if pipeline.cas_col and pipeline.cas_col in all_cols: rec.append(pipeline.cas_col)
            if pipeline.name_col and pipeline.name_col in all_cols: rec.append(pipeline.name_col)

            # Unit and Target columns
            if pipeline.unit_col and pipeline.unit_col in all_cols: rec.append(pipeline.unit_col)
            if pipeline.target_col and pipeline.target_col in all_cols:
                rec.append(pipeline.target_col)

            if st.session_state.dedup_condition_cols:
                for col in st.session_state.dedup_condition_cols:
                    if col in all_cols and col not in rec:
                        rec.append(col)

            keywords = ['Canonical SMILES', 'InChI', 'Is Valid', 'from_web', 'Substructure_', 'Deduplication Strategy']
            for col in all_cols:
                if any(k in col for k in keywords):
                    if col not in rec: rec.append(col)
            return rec

        def select_all():
            st.session_state.export_selected_cols = all_cols

        def select_none():
            st.session_state.export_selected_cols = []

        def select_recommended():
            st.session_state.export_selected_cols = get_recommended_cols()

        c_btn1, c_btn2, c_btn3, c_space = st.columns([1, 1, 1.5, 3])
        c_btn1.button("Select All", on_click=select_all, width="stretch")
        c_btn2.button("Clear All", on_click=select_none, width="stretch")
        c_btn3.button("✨ Recommended", on_click=select_recommended, width="stretch",
                      help="Select key identifiers and processed results only.")

        selected_cols = st.multiselect(
            "Choose columns to include (Drag to reorder):",
            options=all_cols,
            default=st.session_state.export_selected_cols,
            key='export_selected_cols'
        )

        if selected_cols:
            st.caption("Preview:")
            st.dataframe(pipeline.df[selected_cols].head(5), width='content')

            col_gen, col_dl, col_space = st.columns([1, 1, 4])
            msg_container = st.empty()
            with col_gen:
                if st.button("Generate Download File", type="primary", width="stretch"):
                    try:
                        pipeline.save_results(filename, columns=selected_cols)
                        with open(filename, "rb") as f:
                            st.session_state.export_data = f.read()
                            st.session_state.export_fname = filename
                        if os.path.exists(filename): os.remove(filename)
                        msg_container.success(f"✅ File **{filename}** generated successfully! Ready to download ->")
                    except Exception as e:
                        msg_container.error(f"Error: {str(e)}")

            with col_dl:
                file_data = st.session_state.export_data
                file_name = st.session_state.export_fname or filename

                st.download_button(
                    label="⬇️ Click to Download",
                    data=file_data if file_data else b"",
                    file_name=file_name,
                    mime="application/octet-stream",
                    disabled=(file_data is None),
                    type="primary" if file_data else "secondary",
                    width = "stretch"
                )
        else:
            st.warning("⚠️ No columns selected! File will be empty.")

# --- Bottom: Data Preview ---
st.divider()
with st.expander("📜 Processing History (Audit Log)", expanded=False):
    if pipeline.df is not None:
        history_df = pipeline.get_history()
        if not history_df.empty:
            st.dataframe(history_df, width='stretch', hide_index=True)
        else:
            st.caption("No processing steps recorded yet.")
    else:
        st.caption("Load data to start tracking history.")
st.subheader("📊 Data Preview (First 50 Rows)")
if st.session_state.df_preview is not None:
    df_head = st.session_state.df_preview.head(50)
    styled_df = df_head.style

    valid_user_cols = [c for c in st.session_state.user_cols if c in df_head.columns]
    valid_gen_cols = [c for c in st.session_state.gen_cols if c in df_head.columns]

    if valid_user_cols:
        styled_df = styled_df.set_properties(
            subset=valid_user_cols,
            **{'background-color': '#FFF8DC', 'color': 'black'}
        )
    if valid_gen_cols:
        styled_df = styled_df.set_properties(
            subset=valid_gen_cols,
            **{'background-color': '#E6F3FF', 'color': 'black'}
        )
    st.dataframe(styled_df, width="stretch")
else:
    st.caption("No data loaded.")
