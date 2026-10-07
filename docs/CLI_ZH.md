# DiPTox CLI 操作指南

批量合并使用 `--groups '[{"values":["Egg","Embryo"],"replacement":"Embryonic"},{"values":["Fry","Larva"],"replacement":"Larval"}]'`，不再传 `--values`。JSON 流水线使用 `params.groups`。一次生成一列，同时匹配原值；冲突规则报错，未列出的值保留原样。完整示例见 `examples/life_stage_merge_groups.json`。

新增 `merge-values`：选择列和多个值，统一替换为文本标签，生成新列并保留全部行和原列。例如 `--column Species --values '["mouse","rabbit"]' --replacement other --output-column "Species group"`，同时提供 `--input` 和 `--output`。新列可用于后续去重条件，空选择不做修改。JSON 示例见[列调整说明](CONDITION_COLUMNS.md)。

新增 `transform-column`，支持独立指定条件数值/单位列，换算后可选 log10 / -log10。命令及 JSON 示例见[条件列转换](CONDITION_COLUMNS.md)。

[返回项目说明](../README_ZH.md) · [English](CLI.md)

## 按列取值筛选

查看某列的全部不同取值和对应行数：

```bash
python -m diptox column-values --input data.csv --column Species
```

筛选模式 `keep`（仅保留选中值）与 `remove`（去除选中值）互斥。未提供 `--values` 或传入 `[]` 时全保留；JSON `null` 表示缺失值，空字符串 `""` 是独立取值，数值 `1` 与文本 `"1"` 区分。`--values` 接受 JSON 数组，例如数值列：

```bash
python -m diptox filter-values --input data.csv --column Dose --mode remove --values '[0,null]' --output kept.csv --excluded excluded.csv
```

也可使用 JSON 流水线配置字符串选项；多个列筛选步骤按顺序执行，结果为各步保留行的交集。以下内容放入配置的 `steps` 数组：

```json
[
  {"op": "filter-values", "params": {"column": "Species", "mode": "keep", "values": ["rat", "mouse"]}},
  {"op": "filter-values", "params": {"column": "Study type", "mode": "remove", "values": ["in vitro"]}}
]
```

单纯按列筛选不要求 SMILES 列映射。被筛除的整行写入排除记录，包含原因、步骤和源行标识；这是主动筛选，`--strict` 不将其视为数据错误。

Python API 对应 `pipeline.get_column_values(column)` 和 `pipeline.filter_by_values(column, values, mode='keep')`；`pipeline.undo()` 可撤销筛选。GUI 在 **Search & filter → Filter by column values** 中选列、选模式和取值，然后点击 **Apply column filter**；可逐列应用，未选择值则保持数据不变。

## Excel 中的单位显示

CSV/TXT 的负对数单位使用数学减号，例如 `−log10(mol/L)`，并使用 UTF-8 BOM，避免 Excel 将其解释为公式。XLSX 字符串明确保存为文本单元格。重新导入 DiPTox 时，映射的单位列会恢复内部 `-log10(...)` 形式，仍能识别已转换尺度，避免重复取对数。

GUI 的 **Recommended** 下载列会随处理更新，包含当前最终目标值、单位和实际去重条件列。手动修改列选择后保留自定义选择；再次点击 **Recommended** 恢复自动更新。

导出页的 **File name** 可设置下载文件名，扩展名按所选格式自动补全；已输入的常用扩展名会被替换。排除记录使用同一名称加 `-excluded.csv`，例如“筛选结果.xlsx”和“筛选结果-excluded.csv”。命名框只填写文件名，保存目录由浏览器下载设置决定。

除非另有说明，文中的命令均在仓库根目录执行；`examples/cli/...` 指仓库中的示例路径。

CLI 面向智能体和无人值守脚本。安装本版本后使用 `diptox`，也可在仓库内使用 `python -m diptox`；源码安装可执行 `python -m pip install -e .`，随后便可在其他目录调用。`diptox-gui` 仍用于启动图形界面。CLI 支持能力发现、数据检查、预处理、单位转换、去重、子结构搜索、原子数过滤、InChI 计算、化学规则配置和有时限的网络补全，可通过 JSON 流水线按顺序组合处理步骤。

## 检查、配置与执行

在仓库根目录运行：

```bash
python -m diptox schema
python -m diptox inspect --input examples/cli/data.csv --limit 3
python -m diptox run --config examples/cli/pipeline.json --dry-run
python -m diptox run --config examples/cli/pipeline.json
```

[示例配置](../examples/cli/pipeline.json) 处理 [6 条离线记录](../examples/cli/data.csv)，包括两组重复记录、混合单位、一条无效结构和一条无效数值。结果写入 `examples/cli/output/` 下的 `clean.csv`、`excluded.csv` 和 `report.json`。最终保留 2 条记录，标准化数值分别为 1 和 2 mg/L，审计表包含 4 个排除事件。检查产物后，再次运行可加 `--overwrite`。

`schema` 返回配置结构、步骤参数和默认值；`inspect` 返回列名、类型和有限行数的预览。检查时会读取所选数据表，`--limit` 只限制显示行数，不限制读取量。Excel 有多个工作表时，先列出表名，再通过 `--sheet "Sheet1"` 或从 0 开始的 `--sheet-index 0` 选择。

`--dry-run` 检查配置、输入元数据、输入列和生成列、化学规则、查询语法、单位规则是否齐全以及输出路径冲突；网络补全还检查服务能力和所需环境凭据，但不发送网络请求。它会读取输入，但不执行处理步骤、不逐个验证分子、不保证转换或远程查询一定成功，也不写出产物。

## 单步命令与配置

```bash
diptox preprocess --input data.csv --smiles-col SMILES --output processed.csv --excluded rejected.csv --report preprocess.json
diptox units --input data.csv --target-col Value --unit-col Unit --standard-unit mg/L --output units.csv
diptox deduplicate --input processed.csv --smiles-col "Canonical SMILES" --data-type smiles --output unique.csv
```

使用 `diptox 命令 --help` 查看完整参数。列角色包括 `smiles`、`target`、`unit`、`cas`、`name`、`inchikey`、`id`，对应 `--smiles-col`、`--target-col` 等参数；需明确指定步骤要求的列。无表头文件使用 `--no-header`，列名为 `0`、`1` 等，例如 `--smiles-col 0 --target-col 1 --unit-col 2`。文本输入默认 UTF-8；CSV 默认逗号分隔，TXT 默认制表符分隔，可通过 `--delimiter` 指定。

JSON 配置必须包含 `"schema_version": "1"`、`input`、非空且有序的 `steps` 和 `output`。步骤包含 `op` 与 `params`；列映射放在 `input.columns`，产物路径放在 `output.path`、可选的 `excluded_path` 和 `report_path`。未知字段、错误类型和不兼容的设置会在校验时失败。主要默认值如下：

| 设置 | 默认值 / 可选值 |
| --- | --- |
| 输入 | `header: true`、`encoding: "utf-8"` |
| 执行 | `overwrite: false`、`strict: false`；全程不询问交互输入 |
| 预处理 | `mixture_mode: "reject"`、`element_policy: "allow_all"`、`n_jobs: 1`、`chunksize: 100`、`hac_threshold: 3` |
| 化学开关 | 默认启用去盐、去溶剂、去无机物、中和、去同位素、去氢、sanitization 和自由基拒绝；默认关闭去立体、加氢（`add_hs`）、非中性结构拒绝和严格原子检查 |
| 单位转换 | 必须指定 `standard_unit`；分子量来源默认未选；转换需要 MW 时必须选择 `molecular_weight_source` 或 `molecular_weight_col` |
| 去重 | `data_type: "continuous"`、`method: "auto"`、`p_threshold: 0.05`、`log_transform: "None"`、`condition_cols: []`、`dropna_conditions: false` |

布尔步骤参数支持正反形式，例如 `--remove-salts` / `--no-remove-salts`。混合物策略可选 `keep`、`reject`、`largest`；元素策略可选 `allow_all`、`reject_metals`、`allowed_atoms`。连续数据去重通过 `--method auto|IQR|3sigma` 选择异常值筛除方式，通过 `--aggregation mean|max|min` 选择筛除后的聚合方式（JSON 对应 `method` 和 `aggregation`）。默认为 Auto + 平均值；组内不超过 3 条记录时跳过异常值筛除。旧配置中的 `method: "max"` / `"min"` 需将最大/最小值选项移至 `aggregation`；离散数据默认 `vote`，也支持 `priority` 并要求非空优先级列表。仅按结构去重使用 `data_type: "smiles"` 和 `method: "auto"`。对数转换使用字符串 `"None"`、`"-log10"`、`"log10"`，仅适用于连续数据；CLI 中负号开头的值写作 `--log-transform=-log10`。

自定义单位规则格式为 `"conversion_rules": [{"from": "custom", "to": "mg/L", "formula": "x * 1000"}]`；单步命令通过 `--conversion-rules rules.json` 读取同样的数组。公式支持 `x`、可选分子量 `mw`、算术运算和 `log`、`log10`、`exp`。所有使用 `mw` 的转换（包括质量/摩尔浓度转换）都必须明确选择 `molecular_weight_source: "original"` 或 `"standardized"`，也可指定 `molecular_weight_col`。选择 `standardized` 前必须先执行预处理；不再接受分子量来源 `auto`。仅单位倍率换算可不选择分子量来源。

`add_hs: true`（CLI 参数 `--add-hs`）在全部化学处理结束后添加显式氢。可与 `remove_hs: true` 同时开启：先去氢，再在最后加氢。`mixture_mode: "largest"` 遇到并列最大片段时以 `Ambiguous parent` 拒绝，不任意选择片段。Python 接口也建议使用 `mixture_mode` 和 `element_policy`；旧开关仅为兼容保留。

预处理始终去除 `Canonical SMILES` 的原子映射编号，原始列保留来源标记；仅编号不同不计为不同母体结构，也不算化学结构改变。除盐、除溶剂和混合物处理均仅在剩余组分全部相同时去除重复（`A.A → A`），不进行局部去重（`A.A.B → A.B`）。身份比较忽略映射编号，保留该步骤中现存的立体、同位素、电荷和氢表示。全部为同一种已知溶剂的重复项也合并。混合物步骤先合并全部相同的组分，再执行 `keep/reject/largest`，故同一组分的重复项不受最大片段尺寸门槛限制；不同组分仍按所选混合物策略处理。不增加新开关。

乙二醇 `OCCO`、2-甲氧基乙醇 `COCCO` 由 `remove_solvents` 控制，不再进入盐移除规则；作为单独受试分子输入时保留。除盐与除溶剂共用完整片段匹配，兼容常规二碳取代亚砜的 `S=O` 和 `[S+][O-]` 表示。兼容转换仅作用于匹配副本，仍要求整个片段的原子数、键数及查询条件匹配，不把含有该官能团的大分子当成 DMSO 删除，也不执行全局互变异构体或电荷归一化。

`reject_metals` 按碱金属、碱土金属、d 区金属、镧系、锕系及明确列出的 p 区金属判定；该集合与原子白名单独立。贵气体、卤素和 B/Si/Ge/As/Sb/Te 不属于金属集合。边界元素 Po 按本工具约定归为金属；分类及来源见 `diptox/element_policy.py`。金属拒绝理由、盐片段保护及 `Original/Final Metal Elements` 使用同一判定。

配置文件里的相对路径以该配置文件所在目录为基准。`run --config -` 从标准输入读取 JSON，此时路径以当前工作目录为基准。例如先进入 `examples/cli`，在 POSIX shell 中执行 `cat pipeline.json | python -m diptox run --config -`；PowerShell 中执行 `Get-Content -Raw pipeline.json | python -m diptox run --config -`。单步命令的路径也相对当前工作目录。`run --overwrite` 和 `run --strict` 可启用配置中对应的设置。

## 搜索、过滤与标识符

```bash
diptox search --input data.csv --smiles-col SMILES --query-pattern "[OX2H]" --query-type smarts --mode matches --output hydroxyl.csv
diptox search --input data.csv --smiles-col SMILES --query-pattern "CCO" --query-type smiles --mode annotate --output annotated.csv
diptox filter-atoms --input data.csv --smiles-col SMILES --min-heavy-atoms 3 --max-heavy-atoms 30 --max-total-atoms 100 --output filtered.csv
diptox inchi --input data.csv --smiles-col SMILES --output identifiers.csv
```

`search` 必须提供 `query_pattern`；`query_type` 默认为 `smarts`，也可指定 `smiles`。各模式均添加 `Substructure_<query_pattern>` 列：默认的 `annotate` 保留全部行，`matches` 只保留匹配结构，`nonmatches` 只保留不匹配结构。无效或预处理拒绝的结构标记为未知（空值），在两个筛选模式中都会被排除，不计入不匹配结构。

`filter-atoms` 支持 `min_heavy_atoms`、`max_heavy_atoms`、`min_total_atoms`、`max_total_atoms`。边界值是包含端点的非负整数，至少指定一个边界，同类最小值不能大于最大值；总原子数包含氢原子。`inchi` 添加 `InChI` 列，无效结构或计算失败会写入审计。若此前已执行预处理，这些操作使用标准化后的结构。

[离线 P2 示例](../examples/cli/pipeline-p2.json) 复用 `data.csv`，依次执行预处理、含氧子结构筛选、三个重原子的过滤、InChI 计算和仅按结构去重，最终在 `examples/cli/output-p2/` 输出一条乙醇记录，并分别记录无效结构与正常过滤事件。

```bash
python -m diptox run --config examples/cli/pipeline-p2.json --dry-run
python -m diptox run --config examples/cli/pipeline-p2.json
```

## 当前调用的化学规则

```bash
diptox rules
diptox rules --rules-file examples/cli/chemical-rules.json
diptox preprocess --input data.csv --smiles-col SMILES --element-policy allowed_atoms --rules-file examples/cli/chemical-rules.json --output processed.csv
diptox run --config examples/cli/pipeline-p2.json --rules-file examples/cli/chemical-rules.json
```

`rules` 以 JSON 返回有效化学规则。所有处理命令和 `run` 均支持 `--rules-file`：文件直接保存规则变更对象，流水线配置则把同一对象放在顶层 `rules` 字段。命令行规则文件会**整体替换**配置中的 `rules`，不会与其合并。变更只影响当前调用，不修改全局默认值，也不生成持久设置。

```json
{
  "atoms": {"add": ["Xe"], "remove": []},
  "salts": {"add": ["[Na+]"], "remove": []},
  "solvents": {"add": ["CCOCC"], "remove": []},
  "neutralization": {
    "add": [{"reactant": "[O-;X1]", "product": "O"}],
    "remove": []
  }
}
```

各规则组和变更列表均可省略。`atoms` 使用元素符号，`salts` 使用 SMARTS，`solvents` 使用 SMILES；中和规则新增项包含反应物 SMARTS 与替换片段 SMILES，删除项使用反应物 SMARTS 标识规则。先删除再添加；删除不存在的规则会明确报错。预检查会验证化学语法。原子列表变更在预处理使用 `element_policy: "allowed_atoms"` 时生效。报告包含有效规则数量和 SHA-256 指纹，便于复现。

## 网络补全

```bash
diptox enrich --input examples/cli/data-enrich.csv --name-col Name --sources pubchem --send name --request smiles cas iupac mw --timeout 15 --deadline 120 --output enriched.csv --report enrich-report.json
python -m diptox run --config examples/cli/pipeline-enrich.json --dry-run
```

[网络示例](../examples/cli/pipeline-enrich.json) 从名称出发查询结构和属性，再预处理返回的 SMILES 并计算 InChI；去掉 `--dry-run` 才会访问所选服务。产物写入 `examples/cli/output-enrich/`。输入映射只描述原文件：补全步骤请求 `smiles`、`cas` 或 `name` 后，后续步骤可直接使用生成的列角色。如果在预处理后请求 SMILES，则保留既有标准化结构的角色。

| 参数 | 默认值 / 可选值 |
| --- | --- |
| `sources` | `["pubchem"]`；按顺序尝试 `pubchem`、`chemspider`、`comptox`、`cactus`、`chembl`、`cas` |
| `send` | `["smiles"]`；按顺序使用 `smiles`、`cas`、`name`，每个所选角色均需映射列 |
| `request` | 必填非空列表，可选 `smiles`、`cas`、`iupac`、`mw`、`name` |
| `max_workers` | `4`，允许 `1` 至 `32` |
| `timeout` | 每次 HTTP 请求的超时参数，默认 `15` 秒 |
| `deadline` | 补全工作进程总时限，默认 `120` 秒，包含启动、请求、重试与等待 |
| `retries` | 短暂 HTTP 故障的额外尝试次数，默认 `2`；`0` 表示仅发起首次尝试 |
| `retry_delay` | 重试等待间隔，默认 `1` 秒 |
| `interval` | 工作进程内各 HTTP 请求启动的最小间隔，默认 `0.3` 秒 |

CLI 始终直接访问 HTTP API，即使安装了可选的服务 SDK 也不会启用。所选服务不支持某个属性或查询标识类型时，预检查会失败。ChemSpider、CompTox、CAS 的凭据分别从 `DIPTOX_CHEMSPIDER_API_KEY`、`DIPTOX_COMPTOX_API_KEY`、`DIPTOX_CAS_API_KEY` 环境变量读取；选择这些服务前需设置对应变量。配置文件不接受密钥，报告也不写入密钥。例如 PowerShell：

```powershell
$env:DIPTOX_CHEMSPIDER_API_KEY = "your-key"
```

返回属性写入 `<property>_from_web` 列。`Query_Status` 为 `complete`、`partial`、`not_found`、`invalid_identifier` 或 `failed`；`Query_Missing_Fields` 和 `Query_Errors` 保存 JSON 数组；`Data_Source` 保存以属性为键的 JSON 对象，记录服务来源与查询标识类型，`Query_Method` 记录成功使用的标识类型。步骤统计中的 `requested_fields`、`resolved_fields`、`missing_fields` 按“行 × 属性”的单元格计数，并提供 `complete_rows`、`partial_rows`、`failed_rows` 和分类的 `error_counts`。

缺失或部分补全的记录会产生行级警告和审计事件。如果没有任何属性成功返回且发生网络或服务故障，则返回退出码 `4`，不发布产物。超过 `deadline` 会终止工作进程及其请求线程，返回 `DEADLINE_EXCEEDED` 和退出码 `4`，不发布产物；取消操作返回 `130` 并停止工作进程。`--strict` 将不完整的记录结果转为退出码 `5`，在发布前终止；单纯的 `not_found` 不代表网络故障。

## 结果、审计与失败处理

`run` 会在步骤间保留状态和列名变更。本例单位转换生成 `Value (Standardized)`、`Unit (Standardized)`，去重进一步生成 `Value (Standardized)_new`。重新读取 CSV 不会恢复流水线内部状态；分开运行命令时，需显式映射实际生成的列名，包括 `--smiles-col "Canonical SMILES"` 及当前目标值、单位列。可从报告的 `summary.column_roles` 获取最终列角色。

预处理会在主表保留被拒绝的行，并通过 `Is Valid`、`Standardization Status`、`Processing Log` 标记；单位失败也保留并记录状态，去重则移除不可用行并合并重复组。保留列 `Diptox Source Rows` 跟踪从 0 开始的输入记录位置（不含表头），合并后用分号连接。排除表记录每一步的事件，包含 `CLI Step`、`CLI Operation`、`CLI Reason`、`CLI Severity`。同一源记录可能在多个步骤出现，因此 `exclusion_events` 不等于唯一排除记录数。正常的搜索/原子数筛选排除标为 `filtered`，计入 `filtered_rows`，不会导致严格模式失败；无效结构标为 `error`，仍会导致严格模式失败。正常合并的重复行另行统计。

默认情况下，记录级失败产生警告和审计产物，退出码为 0。加 `--strict` 后，此类失败返回退出码 5，且在发布产物前终止，已有文件保持不变。输入、配置和各产物路径必须互不冲突；不指定 `--overwrite` 时禁止覆盖输出。所有文件先写入临时文件，再**逐文件原子发布**，不构成多个文件的整体事务；完成报告最后发布。覆盖时会在发布前移走旧完成报告；若新产物仅部分提交，不会保留旧成功报告。若发布中途失败，前面的文件可能已提交，请结合进程结果和 `error.details.published_artifacts` 判断是否完成。

输入格式：CSV、XLSX、XLS、TXT、SMI、SDF、MOL。主结果格式：CSV、XLSX、TXT、SMI、SDF；不支持 XLS 导出。读取 XLS 需要安装兼容的 pandas Excel 引擎。排除表支持 CSV/XLSX/TXT，报告使用 JSON。Python、GUI 和 CLI 共用 SDF 加载器：失败记录保留为无效行，并附从 1 开始的 `Input Record` 和 `Import Error`。SDF 导出以空分子块保留无效行，标记 `SDF Export Status: Invalid structure placeholder`，原始数据和审计属性仍保留。SMI 导出要求结构有效。可用 `--columns` / `output.columns` 选择输出列，包括预处理生成的 `Original Structure Type`。

SDF 始终以 molblock 结构图为准。非空 SMILES 声明无效或与结构图冲突时，标记 `Structure Reliability=Unreliable`，在导入时移出工作数据，原声明和原因保留在排除表。CLI 将其记录为第 0 步 `load`，严格模式阻止发布。比较忽略原子映射编号和普通显式氢，保留立体、同位素、电荷及组分数量差异。导出从选定结构列重建分子，同步通用 SMILES 属性，并用 `SDF Source SMILES` 归档原始声明。

管线为每行建立稳定的内部唯一身份，避免重复输入索引造成结果或排除记录串行。内部身份不作为界面字段或导出列；展示的来源仍使用原始索引标签，重复去重保留最初的组成员。直接调用低层 `DataDeduplicator` 时需提供唯一 DataFrame 索引；`DiptoxPipeline.load_data` 自动建立内部身份。

命令默认在 stdout 输出一个 UTF-8 JSON 对象，包含 `schema_version`、`status`、`summary`、`artifacts`、`warnings`、`error`。失败时提供 `error.code`、`message`、`details`、`retryable`。日志写入 stderr；`--verbose` 启用信息日志，`--debug` 增加失败堆栈。`--format text` 仅改变响应展示形式，不改变数据文件格式；`--help` 和 `--version` 输出普通文本。

| 退出码 | 含义 |
| --- | --- |
| 0 | 成功，可能包含记录级警告 |
| 1 | 非预期内部错误 |
| 2 | 参数、配置或数据错误 |
| 3 | 文件、产物发布或缺失依赖错误 |
| 4 | 网络/服务故障或补全总时限超时 |
| 5 | 严格模式下记录处理失败，未发布产物 |
| 130 | 操作取消 |

**日志兼容性变化：** 导入 DiPTox 不再自动创建日志目录或修改应用的 root logger。文件日志需显式启用，CLI 默认禁用。原先依赖自动日志文件的 Python 程序可配置 `LogManager`：

```python
from diptox import LogManager

LogManager().configure(enable_file=True)
```
