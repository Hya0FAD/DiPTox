# DiPTox - 计算毒理学的数据整合与清洗

![PyPI Test Version](https://img.shields.io/badge/testpypi-1.0.0-blue) ![License](https://img.shields.io/badge/license-MIT-blue.svg) ![Python Version](https://img.shields.io/badge/python-3.8+-brightgreen.svg)

<p align="center">
  <img src="assets/TOC.png" alt="DiPTox 工作流示意图" width="500">
</p>

**DiPTox** 是一个专为分子数据集的稳健预处理、标准化及多源数据整合而设计的 Python 工具包，专注于计算毒理学工作流。

## 核心功能

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
-   **将分子标准化**为规范的SMILES
-   **按原子数量过滤**（重原子或总原子）

#### 数据去重
-   为重复条目提供灵活的去重策略（`smiles`、`continuous`、`discrete` 目标值）。
-   可自定义的匹配条件和异常值处理方法（`auto`、`IQR`、`3sigma`）。

#### 标识符与属性集成（通过Web服务）
-   从多个在线数据库（**PubChem、ChemSpider、CompTox、Cactus**）查询化学属性。
-   获取并验证 **CAS号、IUPAC名称、InChIKey** 等标识符。
-   通过高性能的**并发请求**加速数据获取。
-   为需要身份验证的服务提供集中的 API 密钥管理。

#### 实用工具
-   使用 SMILES 或 SMARTS 模式执行**亚结构搜索**。
-   **自定义化学处理规则**，包括中和反应、盐/溶剂列表和有效原子。
-   **显示**当前所有生效的处理规则的摘要。

## 安装
```bash
pip install diptox
```

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
  check_valid_atoms=False,    # 检查所有原子是否在有效列表中。默认: False。
  strict_atom_check=False,    # 若为True，则丢弃含无效原子的分子；若为False，则尝试移除它们。默认: False。
  remove_stereo=False,        # 移除立体化学信息 (如 @, / \)。默认: False。
  remove_isotopes=True,       # 移除同位素信息 (如 [13C])。默认: True。
  remove_hs=True              # 移除显式的氢原子。默认: True。
)

# 配置去重
DP.config_deduplicator(condition_cols, data_type, method, custom_method)
DP.data_deduplicate()

# 配置Web查询
DP.config_web_request(source='pubchem/chemspider/comptox/cactus', max_workers, ...)
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

**注意**：`ChemSpider` 和 `CompTox` 需要 API 密钥。请在配置时提供：
```python
DP.config_web_request(
    source='chemspider/comptox',
    chemspider_api_key='your_personal_key',
    comptox_api_key='your_personal_key'
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
- **可选依赖** (根据需要安装):
  - `pubchempy>=1.0.4`: 用于 PubChem 集成
  - `chemspipy>=2.0.0`: 用于 ChemSpider 集成 (需要 API 密钥)
  - `ctx-python>=0.0.1a7`: 用于 CompTox Dashboard 集成 (需要 API 密钥)

## 许可证
本项目采用 MIT 许可证 - 详见 [LICENSE](LICENSE) 文件。

## 支持
如有问题，请在 [GitHub Issues](https://github.com/Hya0FAD/DiPTox/issues) 上提交。