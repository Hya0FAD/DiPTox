# DiPTox - 计算毒理学数据整合与清洗

![PyPI Test Version](https://img.shields.io/badge/testpypi-1.2.0-blue) ![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg) ![Python Version](https://img.shields.io/badge/python-3.8+-brightgreen.svg) [![English](https://img.shields.io/badge/-English-blue.svg)](./README.md)

<p align="center">
  <img src="assets/TOC.png" alt="DiPTox 工作流示意图" width="500">
</p>

**DiPTox** 是一个专为分子数据集的稳健预处理、标准化及多源数据整合而设计的 Python 工具包，专注于计算毒理学工作流。

## v1.2 新特性：交互式 GUI
DiPTox 现在包含一个基于 Streamlit 的用户友好型 **图形用户界面 (GUI)**。用户可以通过可视化的 Web 界面执行数据加载、预处理、Web 检索和去重，无需编写任何代码。

## 核心功能

#### 图形用户界面 (GUI)
-   **可视化操作**：通过浏览器完全控制工作流
-   **实时预览**：应用规则后即时查看数据变化
-   **规则管理**：交互式添加/移除有效原子、盐和溶剂

#### 化学预处理与标准化
一个可配置的管道，用于按特定、可控的顺序清洗和规范化化学结构：
-   **移除盐**
-   **移除溶剂**
-   **处理混合物**（例如，保留最大碎片）
-   **移除无机分子**
-   **中和电荷**
-   **验证原子组成**（基于允许的元素列表）
-   **移除显式氢**
-   **移除立体化学信息**
-   **移除同位素**
-   **移除自由基**
-   **将分子标准化**为规范的SMILES
-   **按原子数量过滤**（重原子或总原子）

#### 数据去重
-   为重复条目提供灵活的去重策略（基于 `smiles` 或基于 `continuous`/`discrete` 目标值）。
-   可自定义的匹配条件（如温度、压强）和去重处理方法（`auto`、`IQR`、`3sigma` 或自定义方法）。

#### 标识符与属性集成（通过Web服务）
-   从多个在线数据库（**PubChem、ChemSpider、CompTox、Cactus、CAS Common Chemistry**）获取并互相转换化学标识符（**CAS号、SMILES、IUPAC名称、常用名、分子量**）。
-   通过高性能的**并发请求**加速数据获取。
-   为需要身份验证的服务提供集中的 API 密钥管理。

#### 实用工具
-   使用 SMILES 或 SMARTS 模式执行**亚结构搜索**。
-   **自定义化学处理规则**，包括中和反应、盐/溶剂列表和有效原子。
-   **显示**当前所有生效的处理规则的摘要。

## 安装
```bash
pip install -i https://test.pypi.org/simple/ diptox
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

# 初始化处理器
DP = DiptoxPipeline()

# 加载数据（可以来自文件路径、列表或DataFrame）
DP.load_data(input_data='file_path/list/dataframe', smiles_col, target_col, cas_col)

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
  remove_salts=True,          # 移除盐片段。默认: True。
  remove_solvents=True,       # 移除溶剂片段。默认: True。
  remove_mixtures=False,      # 基于片段大小处理混合物。默认: False。
  hac_threshold=3,            # 用于移除片段的重原子数阈值。默认: 3。
  keep_largest_fragment=True, # 在混合物中保留最大的片段。默认: True。
  remove_inorganic=True,      # 移除常见的无机分子。默认: True。
  neutralize=True,            # 中和分子上的电荷。默认: True。
  reject_non_neutral=False,   # 仅保留形式电荷为零的分子。默认：False。
  check_valid_atoms=False,    # 检查所有原子是否在有效列表中。默认: False。
  strict_atom_check=False,    # 若为True，则丢弃含无效原子的分子；若为False，则尝试从支链移除它们。默认: False。
  remove_stereo=False,        # 移除立体化学信息 (如 @, / \)。默认: False。
  remove_isotopes=True,       # 移除同位素信息 (如 [13C])。默认: True。
  remove_hs=True              # 移除显式的氢原子。默认: True。
  reject_radical_species=True # 移除含有游离基原子的分子。默认：True。
)

# 配置去重
DP.config_deduplicator(condition_cols, data_type, method, custom_method, priority)
DP.data_deduplicate()

# 配置Web查询
DP.config_web_request(sources=['pubchem/chemspider/comptox/cactus/cas'], max_workers, ...)
DP.web_request(send='cas', request=['smiles', 'iupac'])

# 亚结构搜索
DP.substructure_search(query_pattern, is_smarts=True)

# 保存结果
DP.save_results(output_path='file_path')
```

## 高级配置

### Web服务集成
DiPTox 支持以下化学数据库：
-   `PubChem`: https://pubchem.ncbi.nlm.nih.gov/
-   `ChemSpider`: https://www.chemspider.com/
-   `CompTox`: https://comptox.epa.gov/dashboard/
-   `Cactus`: https://cactus.nci.nih.gov/
-   `CAS`: https://commonchemistry.cas.org/

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
  - `pubchempy>=1.0.4`: 用于 PubChem 集成
  - `chemspipy>=2.0.0`: 用于 ChemSpider 集成 (需要 API 密钥)
  - `ctx-python>=0.0.1a7`: 用于 CompTox Dashboard 集成 (需要 API 密钥)

## 许可证
本项目采用 Apache 2.0 许可证 - 详见 [LICENSE](LICENSE) 文件。

## 支持
如有问题，请在 [GitHub Issues](https://github.com/Hya0FAD/DiPTox/issues) 上提交。
