# diptox/web_ui.py
import streamlit as st
import pandas as pd
import os
from diptox.core import DiptoxPipeline

st.set_page_config(
    page_title="DiPTox GUI",
    page_icon="🧪",
    layout="wide"
)

# --- Session State ---
if 'pipeline' not in st.session_state:
    st.session_state.pipeline = DiptoxPipeline()
if 'df_preview' not in st.session_state:
    st.session_state.df_preview = None
if 'rule_logs' not in st.session_state:
    st.session_state.rule_logs = []
if 'highlight_cols' not in st.session_state:
    st.session_state.highlight_cols = []

pipeline = st.session_state.pipeline


def update_highlights(new_cols):
    """Helper to append new columns to the highlight list without duplicates."""
    if 'highlight_cols' not in st.session_state:
        st.session_state.highlight_cols = []

    if isinstance(new_cols, str):
        new_cols = [new_cols]

    current = set(st.session_state.highlight_cols)
    for col in new_cols:
        if col:
            current.add(col)
    st.session_state.highlight_cols = list(current)


# --- Side bar navigation ---
st.sidebar.title("DiPTox Control Panel")
step = st.sidebar.radio("Go to Step:",
                        ["Data Loading", "Preprocessing", "Web Requests", "Deduplication", "Search & Filter", "Export"])


# --- Auxiliary function ---
def update_preview():
    if pipeline.df is not None:
        st.session_state.df_preview = pipeline.df
    else:
        st.session_state.df_preview = None


def add_rule_log(message):
    st.session_state.rule_logs.append(message)
    st.toast(message, icon="✅")


# ==============================================================================
# Data Loading
# ==============================================================================
if step == "Data Loading":
    st.header("Data Loading")
    st.markdown("Supports: `.csv`, `.xlsx`, `.xls`, `.txt`, `.sdf`, `.smi`")

    uploaded_file = st.file_uploader("Upload Dataset", type=['csv', 'xlsx', 'xls', 'txt', 'smi', 'sdf'])

    if uploaded_file:
        temp_filename = f"temp_upload_{uploaded_file.name}"

        with open(temp_filename, "wb") as f:
            f.write(uploaded_file.getbuffer())

        selected_sheet = None
        if uploaded_file.name.endswith(('.xlsx', '.xls')):
            try:
                with pd.ExcelFile(temp_filename) as xl:
                    if len(xl.sheet_names) > 1:
                        st.info(f"Multi-sheet Excel detected. Please select one.")
                        selected_sheet = st.selectbox("Select Excel Sheet:", xl.sheet_names)
                    else:
                        selected_sheet = xl.sheet_names[0]
            except Exception as e:
                st.error(f"Error reading Excel file structure: {e}")

        with st.expander("Column Configuration", expanded=True):
            c1, c2 = st.columns(2)
            smiles_col = c1.text_input("SMILES Column (Optional)", value="", placeholder="e.g., Smiles")
            target_col = c2.text_input("Target Column (Optional)", value="")

            c3, c4 = st.columns(2)
            cas_col = c3.text_input("CAS Column (Optional)", value="")
            id_col = c4.text_input("ID Column (for .smi)", value="")

            name_col = st.text_input("Name Column (Optional)", value="")

        if st.button("Load Data", type="primary"):
            try:
                kwargs = {}
                if selected_sheet:
                    kwargs['sheet_name'] = selected_sheet

                final_smiles_col = smiles_col.strip() if smiles_col.strip() else None
                final_target_col = target_col.strip() if target_col.strip() else None
                final_cas_col = cas_col.strip() if cas_col.strip() else None
                final_name_col = name_col.strip() if name_col.strip() else None
                final_id_col = id_col.strip() if id_col.strip() else None

                pipeline.load_data(
                    input_data=temp_filename,
                    smiles_col=final_smiles_col,
                    cas_col=final_cas_col,
                    target_col=final_target_col,
                    name_col=final_name_col,
                    id_col=final_id_col,
                    **kwargs
                )

                st.session_state.highlight_cols = [
                    c for c in [final_smiles_col, final_target_col, final_cas_col, final_name_col, final_id_col]
                    if c is not None
                ]

                update_preview()
                st.success(f"Successfully loaded {len(pipeline.df)} records!")

            except Exception as e:
                st.error(f"Error loading data: {str(e)}")
            finally:
                try:
                    import time
                    import gc

                    gc.collect()

                    if os.path.exists(temp_filename):
                        os.remove(temp_filename)
                except PermissionError:
                    print(f"Warning: Could not delete temp file {temp_filename} immediately due to file lock.")
                except Exception as e:
                    print(f"Warning: Error deleting temp file: {e}")

# ==============================================================================
# Preprocessing (Rules & Run)
# ==============================================================================
elif step == "Preprocessing":
    st.header("Preprocessing & Standardization")

    if pipeline.df is None:
        st.warning("Please load data first.")
    else:
        st.markdown("### A. Rule Management (Optional)")

        with st.expander("🛠️ Manage & View Chemical Rules", expanded=False):
            st.caption("Add or remove rules. Changes are applied immediately.")

            # 1. Atoms
            c1, c2, c3 = st.columns([3, 1, 1], vertical_alignment="bottom")
            atom_input = c1.text_input("Valid Atoms", placeholder="e.g., Si, Zr", key="atom_in")
            if c2.button("➕ Add", key="btn_add_atom", width="stretch"):
                if atom_input:
                    pipeline.manage_atom_rules(atoms=[x.strip() for x in atom_input.split(',')], add=True)
                    add_rule_log(f"Added atom(s): {atom_input}")
            if c3.button("➖ Del", key="btn_del_atom", width="stretch"):
                if atom_input:
                    pipeline.manage_atom_rules(atoms=[x.strip() for x in atom_input.split(',')], add=False)
                    add_rule_log(f"Removed atom(s): {atom_input}")

            # 2. Salts
            c1, c2, c3 = st.columns([3, 1, 1], vertical_alignment="bottom")
            salt_input = c1.text_input("Custom Salts (SMARTS)", placeholder="e.g., [Hg+2]", key="salt_in")
            if c2.button("➕ Add", key="btn_add_salt", width="stretch"):
                if salt_input:
                    pipeline.manage_default_salt(salts=[x.strip() for x in salt_input.split(',')], add=True)
                    add_rule_log(f"Added salt: {salt_input}")
            if c3.button("➖ Del", key="btn_del_salt", width="stretch"):
                if salt_input:
                    pipeline.manage_default_salt(salts=[x.strip() for x in salt_input.split(',')], add=False)
                    add_rule_log(f"Removed salt: {salt_input}")

            # 3. Solvents
            c1, c2, c3 = st.columns([3, 1, 1], vertical_alignment="bottom")
            solv_input = c1.text_input("Custom Solvents (SMILES)", placeholder="e.g., CCC", key="solv_in")
            if c2.button("➕ Add", key="btn_add_solv", width="stretch"):
                if solv_input:
                    pipeline.manage_default_solvent(solvents=[x.strip() for x in solv_input.split(',')], add=True)
                    add_rule_log(f"Added solvent: {solv_input}")
            if c3.button("➖ Del", key="btn_del_solv", width="stretch"):
                if solv_input:
                    pipeline.manage_default_solvent(solvents=[x.strip() for x in solv_input.split(',')], add=False)
                    add_rule_log(f"Removed solvent: {solv_input}")

            # 4. Neutralization
            st.markdown("**Neutralization Rules**")
            c1, c2, c3, c4 = st.columns([1.5, 1.5, 0.5, 0.5], vertical_alignment="bottom")
            react = c1.text_input("Reactant (SMARTS)", key="neu_r")
            prod = c2.text_input("Product (SMILES)", key="neu_p")
            if c3.button("➕", key="btn_add_neu", help="Add Rule", width="stretch"):
                if react and prod:
                    pipeline.add_neutralization_rule(react, prod)
                    add_rule_log(f"Added rule: {react} -> {prod}")
            if c4.button("➖", key="btn_del_neu", help="Remove Rule by Reactant", width="stretch"):
                if react:
                    pipeline.remove_neutralization_rule(react)
                    add_rule_log(f"Removed rule: {react}")

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

        col_conf1, col_conf2 = st.columns(2)
        with col_conf1:
            st.caption("Cleaning")
            rm_salts = st.checkbox("Remove Salts", True)
            rm_solvs = st.checkbox("Remove Solvents", True)
            rm_mixt = st.checkbox("Remove Mixtures", False)
            rm_inorg = st.checkbox("Remove Inorganic", True)
            rm_stereo = st.checkbox("Remove Stereochemistry", False)
            rm_hs = st.checkbox("Remove Hydrogens", True)
            rm_iso = st.checkbox("Remove Isotopes", True)
        with col_conf2:
            st.caption("Normalization")
            neu = st.checkbox("Neutralize Charges", True)
            rej_rad = st.checkbox("Reject Radicals", True)
            chk_valid = st.checkbox("Check Valid Atoms", False)
            strict = st.checkbox("Strict Atom Check", False)
            sanitize = st.checkbox("RDKit Sanitize", True)
            keep_large = st.checkbox("Keep Largest Fragment", True)

        st.divider()
        c_sub1, c_sub2 = st.columns(2)
        hac_thres = c_sub1.number_input("HAC Threshold (Mixtures)", 3)
        calc_inchi = c_sub2.toggle("Calculate InChI", True)

        if st.button("🚀 Run Preprocessing", type="primary", width="stretch"):
            progress_bar = st.progress(0.0)
            status_text = st.empty()


            def update_progress(current, total):
                if total > 0:
                    progress_bar.progress(min(current / total, 1.0))
                    status_text.text(f"Processing: {current}/{total}")


            with st.spinner("Processing molecules..."):
                try:
                    pipeline.preprocess(
                        remove_salts=rm_salts, remove_solvents=rm_solvs, remove_mixtures=rm_mixt,
                        remove_inorganic=rm_inorg, neutralize=neu, reject_non_neutral=False,
                        check_valid_atoms=chk_valid, strict_atom_check=strict, remove_stereo=rm_stereo,
                        remove_isotopes=rm_iso, remove_hs=rm_hs, keep_largest_fragment=keep_large,
                        hac_threshold=hac_thres, sanitize=sanitize, reject_radical_species=rej_rad,
                        progress_callback=update_progress
                    )

                    new_cols = ['Canonical SMILES', 'Is Valid', 'Processing Log']

                    if calc_inchi:
                        pipeline.calculate_inchi()
                        new_cols.append('InChI')

                    update_highlights(new_cols)

                    update_preview()
                    st.success("Preprocessing complete!")
                except Exception as e:
                    st.error(f"Preprocessing failed: {str(e)}")

# ==============================================================================
# Web Requests
# ==============================================================================
elif step == "Web Requests":
    st.header("Web Data Enrichment")

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

            send_type = st.multiselect(
                "Send Identifier Input (Priority Order)",
                ["smiles", "cas", "name"],
                default=["smiles"]
            )
            if send_type:
                st.caption(f"📥 Input Priority: {' → '.join(send_type)}")
            else:
                st.caption("⚠️ No input type selected")

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
            req_types = st.multiselect(
                "Request Properties (Output)",
                ["smiles", "cas", "name", "iupac", "mw"],
                default=["cas", "iupac"]
            )

            max_workers = st.slider("Max Workers (Threads)", 1, 16, 4, help="Number of concurrent requests.")

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
            with c_adv6:
                st.markdown("Force API Mode",
                            help="ON: Bypass SDKs (like PubChemPy) and force raw API calls.\nOFF: Use SDKs when available (Recommended).")
                force_api = st.toggle("Force API Mode", value=False, label_visibility="collapsed")

        if st.button("Start Web Request", type="primary"):
            if not sources:
                st.error("Please select at least one data source.")
            elif not send_type:
                st.error("Please select at least one input identifier type.")
            else:
                progress_bar = st.progress(0.0)
                status_text = st.empty()


                def update_progress(current, total):
                    if total > 0:
                        percent = min(current / total, 1.0)
                        progress_bar.progress(percent)
                        status_text.text(f"Querying batch: {current}/{total}")
                with st.spinner(f"Querying... This may take a while. (Check terminal for logs)"):
                    try:
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
                            force_api_mode=force_api
                        )
                        pipeline.web_request(send=send_type, request=req_types, progress_callback=update_progress)
                        new_web_cols = [f"{prop}_from_web" for prop in req_types]
                        new_web_cols.extend(['Query_Status', 'Data_Source', 'Query_Method'])
                        update_highlights(new_web_cols)
                        update_preview()
                        st.success("Web request finished.")
                    except Exception as e:
                        st.error(f"Web request failed: {str(e)}")

# ==============================================================================
# Deduplication
# ==============================================================================
elif step == "Deduplication":
    st.header("Data Deduplication")

    if pipeline.df is None:
        st.warning("Please load data first.")
    else:
        avail_cols = list(pipeline.df.columns)

        condition_cols = st.multiselect("Condition Columns (e.g., pH, Temperature)", avail_cols)
        data_type = st.selectbox("Data Type", ["continuous", "discrete", "smiles"])

        method = "auto"
        priority_list = None

        if data_type == "continuous":
            method = st.selectbox("Deduplication Method", ["auto", "3sigma", "IQR"])
        elif data_type == "discrete":
            method = "vote"
            st.info("Priority Rule: If specified values exist in a duplicate group, they are selected first.")
            priority_input = st.text_area("Priority Values (comma separated, e.g., Active, Intermediate)")
            if priority_input:
                priority_list = [x.strip() for x in priority_input.split(",") if x.strip()]

        if st.button("Run Deduplication", type="primary"):
            progress_bar = st.progress(0.0)
            status_text = st.empty()

            def update_progress(current, total):
                if total > 0:
                    percent = min(current / total, 1.0)
                    progress_bar.progress(percent)
                    status_text.text(f"Processing Groups: {current}/{total}")
            try:
                pipeline.config_deduplicator(
                    condition_cols=condition_cols if condition_cols else None,
                    data_type=data_type,
                    method=method,
                    priority=priority_list
                )
                pipeline.dataset_deduplicate()
                cols_to_light = ['Deduplication Strategy']
                if pipeline.target_col and pipeline.target_col.endswith("_new"):
                    cols_to_light.append(pipeline.target_col)
                update_highlights(cols_to_light)
                update_preview()
                st.success("Deduplication complete.")
            except Exception as e:
                st.error(f"Deduplication failed: {str(e)}")

# ==============================================================================
# Search & Filter
# ==============================================================================
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

            if st.button("Search"):
                try:
                    pipeline.substructure_search(query, is_smarts=is_smarts)
                    update_highlights(f'Substructure_{query}')
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

            if st.button("Apply Filter"):
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

# ==============================================================================
# Export
# ==============================================================================
elif step == "Export":
    st.header("Export Results")

    if pipeline.df is None:
        st.warning("No data to export.")
    else:
        st.markdown(f"**Current Dataset Shape:** {pipeline.df.shape[0]} rows, {pipeline.df.shape[1]} columns")

        col_file1, col_file2 = st.columns(2)
        with col_file1:
            save_fmt = st.selectbox("File Format", ["csv", "xlsx", "txt", "sdf", "smi"])
        with col_file2:
            filename = st.text_input("Output Filename", value=f"diptox_processed.{save_fmt}")

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

            tgt = pipeline.target_col
            if tgt:
                if tgt.endswith("_new") and tgt in all_cols:
                    rec.append(tgt)
                elif tgt in all_cols:
                    rec.append(tgt)

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
            st.caption("Preview of export file (First 5 rows):")
            st.dataframe(pipeline.df[selected_cols].head(5), width="stretch", hide_index=True)
        else:
            st.warning("⚠️ No columns selected! File will be empty.")

        st.markdown("###")  # Spacer
        if st.button("Generate Download File", type="primary", width="stretch", disabled=not selected_cols):
            try:
                pipeline.save_results(filename, columns=selected_cols)

                with open(filename, "rb") as f:
                    file_data = f.read()

                st.download_button(
                    label="⬇️ Click to Download",
                    data=file_data,
                    file_name=filename,
                    mime="application/octet-stream",
                    width="stretch"
                )

                if os.path.exists(filename):
                    os.remove(filename)
            except Exception as e:
                st.error(f"Export failed: {str(e)}")

# --- Bottom: Data Preview ---
st.divider()
st.subheader("📊 Data Preview (First 50 Rows)")
if st.session_state.df_preview is not None:
    df_head = st.session_state.df_preview.head(50)
    valid_highlights = [
        col for col in st.session_state.highlight_cols
        if col in df_head.columns
    ]

    if valid_highlights:
        styled_df = df_head.style.set_properties(
            subset=valid_highlights,
            **{
                'background-color': '#FFF8DC',
                'color': '#333333',
                'border-color': '#ffffff'
            }
        )
        st.dataframe(styled_df, width="stretch")
    else:
        st.dataframe(df_head, width="stretch")
else:
    st.caption("No data loaded.")