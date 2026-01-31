# DiPTox - 计算毒理学数据整合与清洗

[![PyPI](https://img.shields.io/pypi/v/diptox)](https://pypi.org/project/diptox/) ![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg) ![Python Version](https://img.shields.io/badge/python-3.8+-brightgreen.svg) [![English](https://img.shields.io/badge/-English-blue.svg)](./README.md)

<p align="center">
  <img src="assets/TOC.png" alt="DiPTox 工作流示意图" width="500">
</p>

**DiPTox** 是一个专为分子数据集的稳健预处理、标准化及多源数据整合而设计的 Python 工具包，专注于计算毒理学工作流。

## 正式版 v1.0 发布
我们很高兴宣布 DiPTox 在 PyPI 上发布了首个正式稳定版本！这一里程碑带来了生产级的稳定性和显著的性能提升：

* **多进程加速 (Multi-Process Acceleration)**：
    * 通过 `n_jobs` 参数将化学预处理任务的速度提升 **10倍以上**。
    * 针对大规模数据集，智能分配任务至多个 CPU 核心。
* **跨平台健壮性 (Cross-Platform Robustness)**：
    * 针对 Windows 多进程环境实现了专门的 **"卫士机制 (Guard Mechanism)"**，有效防止内存爆炸和递归进程死循环问题。
    * 已在 Windows、Linux 和 macOS 环境下通过稳定性验证。
* **增强的数据加载**：
    * 切换至二进制流解析模式读取 `.sdf` 和 `.mol` 文件，彻底解决编码崩溃问题（如 `utf-8` 与 `latin-1` 混淆）。
    * 自动分子结构解析：即使文件中缺少属性列，也能直接从结构块生成 SMILES。

## DiPTox 社区登记 (可选)
为了更好地了解用户群体并改进软件，DiPTox 在首次使用时会提供一个一次性的、可选的用户信息登记。
-   **完全自愿**：您只需点击一下即可跳过。
-   **注重隐私**：收集的信息仅用于学术影响力评估，绝不会被分享。

## 核心功能

#### 1. 图形用户界面 (GUI)
基于 Streamlit 构建，允许用户通过可视化方式执行所有工作流，无需编写代码。
-   **可视化操作**：通过浏览器完全控制工作流。
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
-   **去重方法**：`auto`（自动）、`IQR`（四分位距）、`3sigma`（标准差）、`vote`（投票）或自定义优先级规则。
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
从 PyPI 安装正式稳定版：
```bash
pip install diptox
```

## 图形用户界面 (GUI)
安装完成后，您可以直接从终端启动图形界面：
```bash
diptox-gui
```
该命令将自动在您的默认 Web 浏览器中打开 DiPTox 界面。

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
      remove_mixtures=False,        # 基于片段大小处理混合物。默认: False。
      hac_threshold=3,              # 用于移除片段的重原子数阈值。默认: 3。
      keep_largest_fragment=True,   # 在混合物中保留最大的片段。默认: True。
      remove_inorganic=True,        # 移除常见的无机分子。默认: True。
      neutralize=True,              # 中和分子上的电荷。默认: True。
      reject_non_neutral=False,     # 仅保留形式电荷为零的分子。默认：False。
      check_valid_atoms=False,      # 检查所有原子是否在有效列表中。默认: False。
      strict_atom_check=False,      # 若为True，则丢弃含无效原子的分子；若为False，则尝试从支链移除它们。默认: False。
      remove_stereo=False,          # 移除立体化学信息 (如 @, / \)。默认: False。
      remove_isotopes=True,         # 移除同位素信息 (如 [13C])。默认: True。
      remove_hs=True,               # 移除显式的氢原子。默认: True。
      reject_radical_species=True,  # 移除含有游离基原子的分子。默认：True。
      n_jobs=4                      # 使用 4 个 CPU 核心加速。默认：1.
    )

    # 配置去重与单位标准化
    conversion_rules = {('g/L', 'mg/L'): 'x * 1000', 
                        ('ug/L', 'mg/L'): 'x / 1000',}
    DP.config_deduplicator(condition_cols, data_type, method, custom_method, priority, standard_unit, conversion_rules, log_transform)
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

## 高级配置

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
  - `streamlit>=1.0.0` (运行 GUI 所需)
- **可选依赖** (根据需要安装，如不安装则使用`requests`发送请求):
  - `pubchempy>=1.0.5`: 用于 PubChem 集成
  - `chemspipy>=2.0.0`: 用于 ChemSpider 集成 (需要 API 密钥)
  - `ctx-python>=0.0.1a10`: 用于 CompTox Dashboard 集成 (需要 API 密钥)

## 许可证
本项目采用 Apache 2.0 许可证 - 详见 [LICENSE](LICENSE) 文件。

## 支持
如有问题，请在 [GitHub Issues](https://github.com/Hya0FAD/DiPTox/issues) 上提交。
