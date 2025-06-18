# 安装
# pip install -i https://test.pypi.org/simple/ molppc
# 要求：python>=3.8, rdkit>=2023.3, requests, tqdm, scipy, openpyxl,
# 非必须库：pubchempy, ctx-python, chemspipy
# 目前最新版：0.11.3

# 调用1：文件处理全流程
from molppc import MolecularProcessor

a = MolecularProcessor()
a.load_data(input_data=r"FileName.xlsx/.xls/.csv/.txt/.sdf/.smi",
            smiles_col='Name', target_col='Name', cas_col='Name', id_col='Name')
# smiles_col, target_col, cas_col, id_col(仅限smi文件)都可以有缺省值，非必须
a.preprocess()
a.config_deduplicator(data_type='discrete', method='vote')
# data_type：discrete(离散), continuous(连续)
# method: auto, vote, 3sigma, IQR
# condition_cols：以列表形式放条件列，如['pH','temperature']
# custom_method：放自定义异常值处理函数
a.data_deduplicate()
a.config_web_request(source='comptox/cactus/pubchem', comptox_api_key='your_key', max_workers=4)
a.web_request(send='smiles', request=['cas', 'iupac'])
a.save_results(r"FileName.xlsx/.csv/.txt/.sdf/.smi")

# 调用2：单独的smiles处理全流程
from molppc import MolecularProcessor

b = MolecularProcessor()
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
    'C1=C(C=C(C(=C1O)O)O)C(=O)OC2=CC(=CC(=C2O)O)C(=O)OC[C@@H]3[C@H]([C@@H]([C@H]([C@@H](O3)OC(=O)C4=CC(=C(C(=C4)OC(=O)C5=CC(=C(C(=C5)O)O)O)O)O)OC(=O)C6=CC(=C(C(=C6)OC(=O)C7=CC(=C(C(=C7)O)O)O)O)O)OC(=O)C8=CC(=C(C(=C8)OC(=O)C9=CC(=C(C(=C9)O)O)O)O)O)OC(=O)C1=CC(=C(C(=C1)OC(=O)C1=CC(=C(C(=C1)O)O)O)O)O',
    'C1=CC(=C2C(=C1)OC(O2)(F)F)C3=CNC=C3C#N',
    'C#C',
    'OC(C(O)C(=O)O)C(=O)O.CCC',
    'CO'
]
b.load_data(input_data=smiles, smiles_col='SMILES')
b.remove_neutralization_rule('[$([N-]C=O)]')  # 删除中和规则
b.add_neutralization_rule('[$([N-]C=O)]', 'N')  # 添加中和规则
b.manage_atom_rules(atoms=['Si', 'Zr'], add=True)  # 添加/删除原子识别
b.manage_default_salt(salts=['II', '[Hg+2]', '[Ba+2]'], add=True)  # 添加/删除脱盐的盐类
b.manage_default_solvent(solvents='CCC', add=True)
b.display_processing_rules()
b.preprocess(remove_solvents=True, neutralize=True, remove_inorganic=True, remove_mixtures=False, check_valid_atoms=False, remove_stereo=False, keep_largest_fragment=False)
# remove_salts: 是否脱盐，默认为True
# remove_solvents: 是否除溶剂，默认为True
# remove_mixtures: 是否除混合物，默认为False
# hac_threshold: 混合物去除的重原子阈值，默认为3
# keep_largest_fragment: 是否保留最大片段，默认为True，在去除阈值后的混合物中保留原子最多的一段
# remove_inorganic: 是否去除无机物，默认为True
# neutralize: 是否中性化，默认为True
# check_valid_atoms: 是否检查有效原子，默认为False
# remove_stereo: 是否去除立体化学，默认为False
# remove_hs: 是否去除显式氢，默认为True
b.filter_by_atom_count(min_total_atoms=4)
b.config_deduplicator()  # 不写配置是SMILES去重
b.data_deduplicate()
b.substructure_search(query_pattern='[C@@H]', is_smarts=True)  # 子结构判断
# b.config_web_request(source='comptox/cactus/pubchem', comptox_api_key='your_key', max_workers=4)
# b.web_request(send='smiles', request=['cas', 'iupac'])
# b.save_results('FileName.xlsx')
print(b.df['Canonical smiles'])
print(len(b.chem_processor.remover.salts))

# 调用3：cas获取smiles后处理
from molppc import MolecularProcessor

c = MolecularProcessor()
c.load_data(input_data="FileName.xlsx/.xls/.csv/.txt/.sdf/.smi",
            cas_col='CAS Number', target_col='Value')
c.config_web_request(source='comptox', comptox_api_key='your_key', max_workers=4)
c.web_request(send='cas', request=['smiles', 'iupac'])
c.preprocess()
c.substructure_search(query_pattern='C(=O)O', is_smarts=True)  # 片段匹配
c.config_deduplicator(data_type='continuous')  # 默认为auto
c.data_deduplicate()
c.save_results('FileName.xlsx/.csv/.txt/.sdf/.smi')

# 调用4：只去重smiles，用于预训练
from molppc import MolecularProcessor

d = MolecularProcessor()
d.load_data(input_data="FileName.xlsx/.xls/.csv/.txt/.sdf/.smi",
            smiles_col='Smiles')
d.preprocess()
d.config_deduplicator()
d.data_deduplicate()
d.save_results('FileName.xlsx/.csv/.txt/.sdf/.smi')

# 调用5：输入sdf, smi文件或保存为sdf, smi文件
from molppc import MolecularProcessor

e = MolecularProcessor()
# a.load_data(input_data=r"E:\课题组\数据\降解\AllPublicnew.sdf",
#             smiles_col='SMILES', target_col='ReadyBiodegradability', cas_col='CASRN')
e.load_data(input_data=r"FileName.xlsx/.xls/.csv/.txt/.sdf/.smi",
            smiles_col='smiles', id_col='zinc_id')
e.preprocess()
e.config_deduplicator(data_type='discrete')
e.data_deduplicate()
e.config_web_request(source='comptox', comptox_api_key='your_key', max_workers=4)
e.web_request(send='smiles', request=['cas'])
e.save_results(r"FileName.xlsx/.csv/.txt/.sdf/.smi")

# 其他：自定义异常值检测逻辑
from molppc import MolecularProcessor
import pandas as pd
import numpy as np


def custom_outlier_filter(values: pd.Series):
    clean = values[(values >= -3.5) & (values <= 3.5)]
    if len(clean) > 1:
        return np.mean(clean), "custom_range"
    if len(clean) == 0:
        return np.mean(values), "all"
    return clean, "custom_range"


e = MolecularProcessor()
e.load_data(input_data="FileName.xlsx/.xls/.csv/.txt/.sdf/.smi",
            smiles_col='Smiles', target_col='Value')
e.preprocess(check_valid_atoms=False)
e.config_deduplicator(custom_method=custom_outlier_filter, data_type='continuous')
e.data_deduplicate()
e.save_results('FileName.xlsx/.csv/.txt/.sdf/.smi')
