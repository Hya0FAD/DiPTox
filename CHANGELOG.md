# DiPTox Changelog / 更新日志

[中文 README](README_ZH.md) · [English README](README.md)

记录各版本的重要功能、兼容性变化及影响使用结果的修正。
Highlights of features, compatibility changes, and corrections that affect results.

## 中文

### 1.1.3

- **计算精度修复**：移除单位换算后按原始有效数字进行的中间舍入，目标列和条件列的换算、log/-log 及后续聚合保留浮点精度。旧结果需从原始数值和所选分子量重新计算。例如 0.05 mg/L、MW 336.231 的 pNOEC 从错误的 7.0 修正为约 6.827667748。
- **批量取值合并**：GUI 新增 Multiple groups (JSON) 模式；API 的 `merge_column_values` 新增 `groups` 参数，CLI 新增 `--groups`，支持多组原值一次归类到一个新列。
- **规则校验与追溯**：同时匹配原始值，避免级联替换；冲突规则在修改前报错，未匹配值和全部行保留，支持一次撤销。
- **使用示例**：新增生命阶段五组分类 JSON 和 `UseExample.py` 的 `use_case_9()`。

### 1.1.2

- **条件列转换**：API、GUI、CLI 可指定数值列和单位列，先换算单位，再可选做 log10 / -log10；生成的新列可作为去重条件，异常行进入排除表，支持撤销。
- **列调整与取值合并**：新增独立列调整页面，位于去重上方，集中取值筛选和合并；API、GUI、CLI 可把多个值合并为 `other` 等标签，保留原列和全部行，新列可用于去重，支持撤销。
- **界面布局**：Columns to transform 使用两行等宽布局，窄屏自动改为单列。

### 1.1.1

- **按列取值筛选**：API、GUI 和 CLI 支持查看列内全部取值及计数，按互斥的保留/去除模式筛选；未选值全保留，筛除行进入排除记录，支持撤销。
- **导出修复**：GUI 推荐列随处理更新，包含最终目标值、单位及去重条件；Excel 文本单元格和 CSV/TXT 数学减号避免负对数单位被识别为公式。
- **下载命名**：GUI 导出页可自定义文件名，扩展名随导出格式调整，排除记录使用相同名称加 `-excluded.csv`。

### 1.1.0

- **Python 3.8 兼容**：Python 3.8/3.9 自动安装 NiceGUI 2.24.2，Python 3.10 及以上继续使用 NiceGUI 3.16+；共用新版界面和处理功能，适配文件上传、后台线程、类型注解、命令行布尔参数及 pandas 2.0 的可空数值换算。
- **新版图形界面**：采用 NiceGUI 本地浏览器界面，切换页面时保留配置；耗时操作在后台运行，成功后更新结果。
- **批量处理**：提供命令行流水线、配置预检和运行审计，支持通过 `n_jobs` 并行预处理分子。
- **网络补全与规则配置**：支持多源查询、限流重试、总时限及字段来源记录，可通过 JSON 配置每次运行的化学规则。
- **分子预处理**：支持混合物保留、剔除或唯一最大组分选择，以及元素策略和显式加氢；Canonical SMILES 去除原子映射编号，保留原始结构供追溯。
- **单位与目标值处理**：使用分子量的换算须明确选择原始结构、标准化结构或指定 MW 列，取消 MW 来源的 `auto` 模式；对数变换保留原值和目标尺度，避免重复变换。去重方法中的 `auto` 不变。
- **SDF 结构一致性**：以文件中的分子结构图为准，将与声明 SMILES 冲突的记录移入排除表；导出时同步结构字段并保留原始声明。

### 1.0.6

* **数据读取修复**：修复并优化了对原生 `.smi` (SMILES) 文件的读取逻辑，解决了先前版本中的解析问题，确保大规模化学数据库的稳定加载。
* **Web Request 网络模块全面重构**：针对大规模并发请求进行了加固。引入“能力清单”与针对语义/鉴权错误的短路拦截逻辑（Fast-Fail），消灭了无效重试死循环；同时清除了数据抓取中的“静默失败”现象，将极高颗粒度的精确死因（如 `Failed -> pubchem: Not Found | chemspider: Auth Error (401)`）以及精确到“字段级”的数据溯源记录直接写入结果中，使 API 调试完全透明，极大提升了毒理学数据集的审计置信度。

### 1.0.5

* **增强的单位标准化**：新增对常用数学符号 `^`（幂运算）的支持并自动映射为 `**`， 同时修复了当数据集中仅存在单一单位时系统会强制跳过转换程序的逻辑漏洞。
* **去重逻辑功能升级**：在原有的 `-log10` 基础上新增了 `log10` 转换选项， 使工具包不仅能处理毒性数据（pIC50），还能完美适配水溶解度（logS）或分配系数等理化性质的去重需求。
* **系统健壮性与容错处理**：在单位标准化与去重模块中引入了强制数值校验， 能够自动剔除目标列中的非法字符串（如 "N/A" 或 ">100"）并弹出警告，显著提升了处理真实实验数据的稳定性。
* **关键状态重置修复**：修复了 `load_data` 方法未重置预处理标志位的问题， 确保了用户在重新加载数据集后，网络请求（Web Request）的自动列映射逻辑（如自动识别 `smiles_from_web`）能够恢复正常工作。

### 1.0.4

* **GUI 状态管理修复**：修复了在导出页面点击“撤销上一步”时触发的 `StreamlitAPIException`报错。通过引入 `on_click` 回调函数，在 UI 重新渲染前安全地更新组件状态，确保撤销操作不会崩溃。
* **预处理规则优化**：调整并优化了部分默认的电荷中和规则（`Neutralization rules`）。

### 1.0.3

* **增强的单位标准化**：自定义转换公式现在完全支持分子量（`mw`），可以在摩尔浓度和质量浓度之间进行转换（例如，使用 `x * mw * 1000` 这样的公式）。
* **GUI 界面优化**：Streamlit 图形界面经过了重新设计，呈现出更整洁的布局，减少了视觉干扰，直观地对配置面板进行了分组，并改善了组件对齐。
* **全面的审计记录 (History)**：处理历史记录得到了大幅升级。它现在会详细记录每个操作的颗粒化参数——包括具体触发了哪些预处理规则、激活的去重条件、网络查询状态以及子结构搜索的匹配数量。

## English

### 1.1.3

- **Numeric precision fix**: Removed intermediate rounding to source significant figures after unit conversion. Target/condition conversions, log transforms and downstream aggregation retain floating-point precision. Existing results must be recalculated from original values and the chosen molecular-weight basis. For 0.05 mg/L and MW 336.231, pNOEC is approximately 6.827667748 rather than 7.0.
- **Batch value merging**: Added GUI Multiple groups (JSON) mode, the API `merge_column_values(groups=...)` parameter, and CLI `--groups` to apply multiple category mappings into one new column.
- **Validation and traceability**: Rules match original values simultaneously without cascading. Conflicts are rejected before mutation; unmatched values and all rows are preserved, with single-step undo.
- **Examples**: Added a five-group life-stage JSON example and `use_case_9()` in `UseExample.py`.

### 1.1.2

- **Condition transformations**: API, GUI and CLI accept independent value/unit columns, convert units before optional log10/-log10, and expose generated columns for grouping. Invalid rows are audited; undo restores the operation.
- **Column adjustments and merging**: A dedicated page above Deduplication groups filtering and merging. API, GUI and CLI merge selected values into a label such as `other`, preserving source columns and all rows. Generated columns support grouping and undo.
- **Interface layout**: Columns to transform uses two rows of equally sized fields and switches to one column on narrow screens.

### 1.1.1

- **Column-value filtering**: API, GUI and CLI list distinct values/counts and support exclusive keep/remove modes. Empty selection retains all rows; removed rows are audited and filters can be undone.
- **Export fixes**: Recommended GUI columns follow final targets, units and deduplication conditions. Explicit XLSX text cells and mathematical minus signs in CSV/TXT prevent negative-log units from being interpreted as formulas.
- **Download names**: The GUI export page accepts a custom file name, supplies the selected format's extension, and names exclusions with the same base plus `-excluded.csv`.

### 1.1.0

- **Python 3.8 compatibility**: Python 3.8/3.9 automatically installs NiceGUI 2.24.2; Python 3.10+ continues to use NiceGUI 3.16+. Both share the current interface and processing features, with adapters for uploads, background threads, type annotations, CLI boolean flags, and nullable numeric conversions on pandas 2.0.
- **New graphical interface**: A local NiceGUI browser interface preserves settings across pages. Long-running operations run in the background and update results when successful.
- **Batch processing**: Command-line pipelines support configuration preflight and run reports, with parallel molecular preprocessing through `n_jobs`.
- **Network enrichment and rule configuration**: Supports multiple sources, rate limits, retries, total deadlines, and field provenance. Chemical rules can be configured for each run through JSON.
- **Molecular preprocessing**: Supports retaining or rejecting mixtures, selecting a unique largest component, element policies, and explicit hydrogen addition. Canonical SMILES omits atom-map numbers while retaining the original structure for reference.
- **Units and target values**: MW-dependent conversions require an explicit choice of original structure, standardized structure, or a supplied MW column; the MW-source `auto` mode is removed. Log transformations preserve original values and track their scale to avoid repeated transformation. The deduplication method named `auto` is unchanged.
- **SDF structure consistency**: Uses the molecular graph as the source of truth and moves records with conflicting declared SMILES to the exclusions table. Export synchronizes structure fields and preserves original declarations.

### 1.0.6

* **Data Loading Fixes**: Fixed and optimized the native parsing logic for `.smi` (SMILES) files, resolving previous reading issues to ensure stable ingestion of large-scale chemical databases.
* **Web Request Module Overhaul**: Completely refactored the network request engine for stability and transparency. This update introduces a "Capability Map" and fast-fail logic to intelligently intercept unsupported queries and Auth/404 errors (eliminating infinite retry deadlocks). Furthermore, it eradicates "silent failures" by logging highly granular failure reasons (e.g., `Failed -> pubchem: Not Found | chemspider: Auth Error (401)`), and implements field-level data provenance to strictly record the exact source for each molecular property, drastically improving dataset auditability.

### 1.0.5

* **Enhanced Unit Standardization**: Added support for the standard math operator `^` (power) by automatically mapping it to `**`, and fixed a logic error that caused single-unit datasets to be skipped even when a different target unit was specified.
* **Deduplication Logic Upgrades**: Introduced a `log10` transformation mode alongside the existing `-log10` option, enabling support for both toxicity data (pIC50) and physicochemical properties like water solubility (logS) or partition coefficients.
* **Robustness & Error Handling**: Implemented strict numerical validation using `errors='coerce'` in standardization and deduplication modules to automatically filter out invalid strings (e.g., "N/A", ">100") with clear warning feedback in the GUI.
* **Critical State Management Fix**: Resolved an issue where `load_data` failed to reset the `_preprocess_key` flag, ensuring that automatic column mapping logic for Web Requests (like auto-detecting `smiles_from_web`) functions correctly after a new dataset is loaded.

### 1.0.4

* **GUI State Management Fix**: Resolved a `StreamlitAPIException` on the Export page that occurred when using the "Undo Last Step" feature. Implemented proper `on_click` callbacks to safely mutate the `session_state` (specifically for `export_selected_cols`) before the UI re-renders, ensuring a crash-free and seamless undo experience.
* **Refined Preprocessing Rules**: Adjusted and optimized several default charge neutralization rules.

### 1.0.3

* **Enhanced Unit Standardization**: Custom conversion formulas now fully support molecular weight (`mw`). You can seamlessly convert between molarity and mass concentrations (e.g., using formulas like `x * mw * 1000`).
* **GUI Interface Optimization**: The Streamlit graphical interface has been beautifully redesigned for a more professional, clean, and logical scientific layout. We've reduced visual clutter, grouped configuration panels intuitively, and improved component alignment.
* **Comprehensive Audit Log (History)**: The processing history has been heavily upgraded. It now records granular parameters for every operation—including exactly which preprocessing rules were triggered, active deduplication conditions, web query statuses, and substructure search match counts.
