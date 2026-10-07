# DiPTox - 计算毒理学数据整合与清洗

[![PyPI](https://img.shields.io/pypi/v/diptox)](https://pypi.org/project/diptox/) [![Conda](https://img.shields.io/conda/vn/conda-forge/diptox.svg)](https://anaconda.org/conda-forge/diptox) [![Conda Platforms](https://img.shields.io/conda/pn/conda-forge/diptox.svg)](https://anaconda.org/conda-forge/diptox) ![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg) ![Python Version](https://img.shields.io/badge/python-3.10+-brightgreen.svg) [![English](https://img.shields.io/badge/-English-blue.svg)](./README.md)  [![PyPI下载量](https://static.pepy.tech/badge/diptox)](https://pepy.tech/project/diptox) [![Conda下载量](https://img.shields.io/conda/dn/conda-forge/diptox.svg)](https://anaconda.org/conda-forge/diptox)

<p align="center">
  <img src="assets/TOC.png" alt="DiPTox 工作流示意图" width="500">
</p>

**DiPTox** 是一个专为分子数据集的稳健预处理、标准化及多源数据整合而设计的 Python 工具包，专注于计算毒理学工作流。

## v1.1.3 更新

- **计算精度修复**：单位换算不再按原始有效数字中间舍入，避免 pNOEC 等对数标签失真；旧标签需从原始浓度和分子量重新计算。
- **批量取值合并**：GUI、API、CLI 支持通过 JSON 数组一次配置多组合并规则，只生成一个分组列。所有规则同时匹配原始值，未列出的值及原列保留，冲突规则报错，支持一次撤销。参见[使用说明](docs/CONDITION_COLUMNS.md)及[生命阶段分类示例](examples/life_stage_merge_groups.json)。

## v1.1.2 更新

- **条件列转换**：指定数值列和单位列，换算后可选 log10 / -log10，生成结果可用于去重条件；[API、GUI、CLI 使用说明](docs/CONDITION_COLUMNS.md)。
- **列调整页面**：位于去重上方，集中按列取值筛选和取值合并。可把多个类别合并为 `other` 等标签，生成分组新列并保留原值；API 和 CLI 同步支持。

## v1.1.1 更新

- **按列取值筛选**：API、GUI 和 CLI 可查看取值及计数，选择保留或去除；未选值全保留，筛除行进入排除记录，支持撤销。
- **导出修复**：GUI 推荐列自动跟随最终数值、单位及去重条件；对数单位在 Excel 中按文本显示。
- **下载命名**：GUI 可自定义文件名，扩展名随所选格式调整。

## v1.1.0 更新

- **CLI 与 JSON 流水线**：新增面向智能体和脚本的命令行入口，支持能力发现、数据检查、配置预检，以及预处理、单位转换、去重、子结构搜索、原子数过滤和 InChI 计算。
- **网络补全与化学规则**：支持多源属性查询、限流、重试、总时限与字段来源记录；可通过 JSON 配置当前调用的化学规则。
- **NiceGUI 图形界面**：页面间保留配置，耗时操作在后台运行；修改数据的任务成功后再提交结果。
- **结果与审计**：CLI 提供结构化 JSON 响应、明确的退出码、源行追踪、排除记录和运行报告，便于复现与检查。

操作方法见 [CLI 操作指南](docs/CLI_ZH.md)，当前及以往版本说明见 [更新日志](CHANGELOG.md#中文)。

## DiPTox 社区登记 (可选)
为了更好地了解用户群体并改进软件，DiPTox 在首次使用时会提供一个一次性的、可选的用户信息登记。
-   **完全自愿**：您只需点击一下即可跳过。
-   **注重隐私**：收集的信息仅用于学术影响力评估，绝不会被分享。

## 核心功能

#### 1. 图形用户界面 (GUI)
基于 NiceGUI 构建的本地 Web 界面允许用户通过浏览器执行所有工作流，无需编写代码。
-   **页面状态持久**：在不同处理步骤间切换时，配置值不会丢失。
-   **后台任务响应**：长时间计算和网络请求不会阻塞页面切换。
-   **实时预览**：应用规则后即时查看数据变化。
-   **规则管理**：交互式添加/移除有效原子、盐、溶剂及单位转换公式。
-   **智能列映射**：智能识别表头及二进制文件结构。

#### 2. 化学预处理与标准化
一个可配置的管道，用于清洗和规范化化学结构。
-   **严格的无机物过滤**：更新了 SMARTS 匹配模式，能准确识别复杂的无机物（如离子氰化物）而不误伤有机腈类。
-   **处理流程**：
    -   移除盐与溶剂
    -   处理混合物（保留最大片段）
    -   移除无机分子
    -   电荷中和 & 原子组成验证
    -   移除显式氢、立体化学及同位素信息
    -   **移除自由基**：自动丢弃含有游离基原子的分子。
    -   标准化为规范 SMILES
    -   按原子数过滤

#### 3. 单位标准化
轻松将异构的目标值数据归一化为统一单位。
-   **自动转换**：内置 **浓度**、**时间**、**压力** 和 **温度** 的常用转换规则。
-   **自定义公式**：支持通过 GUI 或脚本交互式定义数学规则（例如 `x * 1000` 或 `10**(-x)`）。
-   **统一输出**：将多样化的单位（如 `ug/mL`, `g/L`, `M`）标准化为单一目标（如 `mg/L`）。

#### 4. 数据去重
提供灵活的重复条目处理策略及高级控制。
-   **数据类型**：支持 `continuous`（连续值，如 IC50）和 `discrete`（离散值，如 Active/Inactive）。
-   **连续值**：`method="auto"`、`"IQR"`、`"3sigma"` 选择异常值筛除方式；`aggregation="mean"`、`"max"`、`"min"` 选择剩余值的平均值、最大值或最小值。默认 Auto + 平均值；组内不超过 3 条记录时跳过异常值筛除。
-   **离散值**：`method="vote"`（投票）或 `"priority"`（优先级）。
-   **Log 变换**：支持在去重逻辑执行**前**应用 `-log10` 变换（例如 IC50 $\to$ pIC50），以正确处理生物活性数据。
-   **灵活的 NaN 处理**：新增选项允许保留条件列中存在缺失值的行（将 *NaN* 视为一个独立分组），防止数据意外丢失。

#### 5. 全流程历史追踪 (审计日志)
-   自动记录 **审计日志 (Audit Log)** 中的每一步操作（加载、预处理、过滤等）。
-   详细追踪 **时间戳**、**操作详情** 以及行数变化（**Delta**）。
-   可通过 API (`get_history()`) 获取或在 GUI 中可视化查看。

#### 6. 标识符与属性集成
-   从多个在线源（**PubChem、ChemSpider、CompTox、Cactus、CAS Common Chemistry、ChEMBL**）获取并互转标识符（**CAS, SMILES, IUPAC, MW**）。
-   具备自动限流与重试机制的高性能 **并发请求**。

#### 7. 实用工具
-   使用 SMILES 或 SMARTS 模式执行**亚结构搜索**。
-   **自定义化学处理规则**，包括中和反应、盐/溶剂列表和有效原子。
-   **显示**当前所有生效的处理规则的摘要。

## 安装
可以通过 `pip` 或 `conda`/`mamba` 安装稳定版：
### 方式1：从 PyPI 安装
从 PyPI 安装官方稳定版本：
```bash
pip install diptox
```
### 方式2：从 Conda-forge 安装
通过以下步骤可以使用 `conda-forge` 频道安装 `diptox`：先将 `conda-forge` 添加到您的软件源列表中。

```bash
conda config --add channels conda-forge
conda config --set channel_priority strict
```

一旦 `conda-forge` 渠道被启用，就可以使用 `conda` 来安装 `diptox`：

```bash
conda install diptox
```

或者通过 `mamba` 安装:

```bash
mamba install diptox
```

## 图形用户界面 (GUI)
安装完成后，您可以直接从终端启动图形界面：
```bash
diptox-gui
```
该命令会在本地启动 DiPTox，并自动在默认浏览器中打开界面。

## 命令行模式（CLI）

CLI 面向智能体和无人值守脚本，支持通过单步命令或 JSON 流水线处理数据。

```bash
diptox --help
```

安装与运行方式、完整参数、示例、网络补全和审计说明见 [CLI 操作指南](docs/CLI_ZH.md)。

## 快速入门
```python
from diptox import DiptoxPipeline

def main():
    # 初始化处理器
    DP = DiptoxPipeline()

    # 加载数据（可以来自文件路径、列表或DataFrame）
    DP.load_data(input_data='file_path/list/dataframe', smiles_col, target_col, cas_col, unit_col)

    # 自定义处理规则（可选）
    print("--- 默认规则 ---")
    DP.display_processing_rules()

    DP.manage_atom_rules(atoms=['Si'], add=True)         # 将 'Si' 添加到有效原子列表
    DP.manage_default_salt(salts=['[Na+]'], add=False)   # 示例：从盐列表中移除钠盐
    DP.manage_default_solvent(solvents='Cl', add=False)  # 示例：从溶剂列表中移除氯
    DP.add_neutralization_rule('[$([N-]C=O)]', 'N')      # 添加一条自定义中和规则

    print("\n--- 自定义后的规则 ---")
    DP.display_processing_rules()

    # 配置预处理流程
    DP.preprocess(
      remove_salts=True,            # 移除盐片段。默认: True。
      remove_solvents=True,         # 移除溶剂片段。默认: True。
      mixture_mode="reject",       # keep / reject / largest。默认: reject。
      hac_threshold=3,              # largest 模式要求唯一最大片段，且重原子数 > 3。
      remove_inorganic=True,        # 移除常见的无机分子。默认: True。
      neutralize=True,              # 中和分子上的电荷。默认: True。
      reject_non_neutral=False,     # 仅保留形式电荷为零的分子。默认：False。
      element_policy="allow_all",  # allow_all / reject_metals / allowed_atoms。默认: allow_all。
                                   # allowed_atoms：含非允许元素即剔除整个分子，不检查原子度数。
      remove_stereo=False,          # 移除立体化学信息 (如 @, / \)。默认: False。
      remove_isotopes=True,         # 移除同位素信息 (如 [13C])。默认: True。
      remove_hs=True,               # 移除显式的氢原子。默认: True。
      add_hs=False,                 # 在全部化学处理结束后添加显式氢。默认: False。
      reject_radical_species=True,  # 移除含有游离基原子的分子。默认：True。
      n_jobs=4                      # 使用 4 个 CPU 核心加速。默认：1.
    )

    # 配置去重与单位标准化
    conversion_rules = {('g/L', 'mg/L'): 'x * 1000', 
                        ('ug/L', 'mg/L'): 'x / 1000',}
    DP.standardize_units(standard_unit="mg/L", conversion_rules=conversion_rules)
    DP.config_deduplicator(condition_cols=condition_cols, data_type=data_type, method=method, aggregation="mean")
    DP.dataset_deduplicate()

    # 配置Web查询
    DP.config_web_request(sources=['pubchem/chemspider/comptox/cactus/cas'], max_workers, ...)
    DP.web_request(send='cas', request=['smiles', 'iupac'])

    # 亚结构搜索
    DP.substructure_search(query_pattern, is_smarts=True)

    # 保存结果
    DP.save_results(output_path='file_path')

    # 查看处理历史
    print(DP.get_history())
    # Output Example:
    #               Step Timestamp  Rows Before  Rows After   Delta                               Details
    # 0     Data Loading  10:00:01            0        1000   +1000                   Source: dataset.csv
    # 1    Preprocessing  10:00:05         1000         950     -50  Valid: 950, Invalid: 50. Order: ...
    # 2    Deduplication  10:00:08          950         800    -150       Method: auto (Log10 Transformed)

# 关键：在 Windows 下使用多进程 (n_jobs > 1) 必须包含此保护块！
# 它可以防止无限递归循环和内存爆炸。
if __name__ == '__main__':
    main()
```

凡公式使用 `mw`，必须明确选择 `molecular_weight_source="original"`、`"standardized"` 或指定 `molecular_weight_col`。`standardized` 要求先完成预处理；分子量来源默认未选择。不使用 MW 的单位倍率换算无需选择。`mixture_mode="largest"` 遇到并列最大片段时以 `Ambiguous parent` 拒绝。当 `remove_hs` 与 `add_hs` 同时开启时，先去氢，并在其他化学处理全部结束后加氢。Python 接口仅为兼容保留旧混合物/元素开关，新代码请使用 `mixture_mode` 和 `element_policy`。

`Canonical SMILES` 去除原子映射编号，原始 SMILES 和 `Original Canonical SMILES` 保留来源标记。除盐、除溶剂和混合物处理均在剩余组分全部相同时只保留一份：`A.A → A`，`A.A.B` 不会局部去重为 `A.B`。混合物模式（包括 `keep`）在此规则之后执行，不增加独立开关。全部为同一种已知溶剂的重复项也会合为一份；其他已选化学规则仍照常执行。

乙二醇 `OCCO` 和 2-甲氧基乙醇 `COCCO` 按溶剂处理，不由盐规则移除。盐和溶剂使用相同的完整片段匹配入口，支持亚砜（例如 DMSO）的 `S=O` / `[S+][O-]` 两种表示。`reject_metals` 使用明确的金属集合；贵气体及 B、Si、Ge、As、Sb、Te 不计为金属，`allowed_atoms` 则另按用户配置的原子列表检查。

## 高级配置

### 连续处理与原始数据保护

输入字段与程序输出列重名时，程序自动将原字段保留为 `原列名 (Input)` 并更新映射，无需预先修改 Excel。SDF 导出也会同步用户指定的自定义结构字段，并归档原始声明。

中和位点按规范顺序选择，避免 SMILES 写法影响结果。启用组分处理时，中和或加氢结束后再次合并全部相同的组分；`A.A.B` 仍不进行局部去重。非字符串 SMILES 按行标记无效，不中断其他记录。

重新预处理改变标准分子量时，依赖旧标准分子量的目标值会清空并标记 `Stale standardized molecular weight`。请从原始浓度重新换算后再去重；已聚合的数据需要先撤销或恢复原始记录。纯倍率换算和不依赖标准结构的 MW 来源不受影响。

`log_transform` 保留原值，将对数结果写入新目标列，单位明确显示为 `-log10(M)` 等形式。相同设置重复去重不会再次取对数。对数目标不能直接用于线性单位换算；先撤销对数处理或重新选择原始线性数据。新指定的去重单位会执行；连续目标中的 `inf/-inf` 会被排除。

### Web服务集成
DiPTox 支持以下化学数据库：
-   `PubChem`: https://pubchem.ncbi.nlm.nih.gov/
-   `ChemSpider`: https://www.chemspider.com/
-   `CompTox`: https://comptox.epa.gov/dashboard/
-   `Cactus`: https://cactus.nci.nih.gov/
-   `CAS`: https://commonchemistry.cas.org/
-   `ChEMBL`: https://www.ebi.ac.uk/chembl/

**注意**：`ChemSpider` 、 `CompTox` 和 `CAS` 需要 API 密钥。请在配置时提供：
```python
DP.config_web_request(
    source='chemspider/comptox',
    chemspider_api_key='your_personal_key',
    comptox_api_key='your_personal_key'
    cas_api_key='your_personal_key'
)
```
## 环境要求
- `Python>=3.8`
- **核心依赖**:
  - `requests`
  - `rdkit>=2023.3`
  - `tqdm`
  - `openpyxl`
  - `scipy`
  - Python 3.8/3.9：`nicegui==2.24.2`；Python 3.10 及以上：`nicegui>=3.16,<4`（安装时自动选择）
- **可选依赖** (根据需要安装，如不安装则使用`requests`发送请求):
  - `pubchempy>=1.0.5`: 用于 PubChem 集成
  - `chemspipy>=2.0.0`: 用于 ChemSpider 集成 (需要 API 密钥)
  - `ctx-python>=0.0.1a10`: 用于 CompTox Dashboard 集成 (需要 API 密钥)

两个 Python 版本分支共用同一套新版 GUI 和处理功能，包括后台任务、化学规则配置、单位换算、去重、筛选与导出。兼容层适配旧版框架的文件上传和 Python 3.8 的后台线程任务。在所需 Python 环境中运行 `python -m pip install .` 即可安装当前源码。Python 3.8/3.9 使用旧版框架，GUI 默认仅监听本机 `127.0.0.1`。

## 许可证
本项目采用 Apache 2.0 许可证 - 详见 [LICENSE](LICENSE) 文件。

## 支持
如有问题，请在 [GitHub Issues](https://github.com/Hya0FAD/DiPTox/issues) 上提交。
