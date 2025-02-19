# 调用1：文件处理全流程
from molppc import MolecularProcessor

a = MolecularProcessor()
a.load_data(input_data=r"E:\周总结\2025\第7次 2025年2月15日\ICICIE_CAR1697.csv",
            smiles_col='SMILES', target_col='Carcinogenicity')
a.process(neutralize=False, remove_stereo=False)
# neutralize: 是否中性化，默认为True
# remove_salts: 是否脱盐，默认为True
# check_valid_atoms: 是否检查有效原子，默认为True
# remove_stereo: 是否去除立体化学，默认为True
# remove_hs: 是否去除隐式氢，默认为True
# keep_largest_fragment: 是否保留最大片段，默认为False
# hac_threshold: 盐去除的重原子阈值，默认为3
# sanitize: 是否执行化学校验，默认为True
a.config_deduplicator(data_type='discrete')  # 离散：discrete，连续：continuous，可选参数method为auto, vote, 3sigma, IQR
a.deduplicate_data()
# a.config_web_request(source='chemspider', chemspider_api_key='key')
# a.config_web_request(source='cactus', max_workers=2)
a.config_web_request(source='pubchem', max_workers=4)
# a.config_web_request(source='comptox', comptox_api_key='key', max_workers=4)
a.add_web_request(cas=True, iupac=True)
a.save_results(r"E:\周总结\2025\第7次 2025年2月15日\ICICIE_CAR1697_.csv")

# 调用2：单独的smiles处理全流程
from molppc import MolecularProcessor

b = MolecularProcessor()
smiles = [
    "C(C(=O)O)C.[Na+]",
    "C1=CC=CC=C1.CCCC",
    "CC1=CC=CC=C1.CCC",
    'CCC(=O)O',
    'CCCC',
    'CC',
    'CCC.CCCC',
    '234',
    'C(Cl)(Cl)(Cl)Cl',
    'Error',
    '[Si](C)(C)[Si](C)(C)C',
    '[N-]C(=O)C',
    'OC[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)CO'
]
b.load_data(input_data=smiles, smiles_col='SMILES')
b.remove_neutralization_rule('[$([N-]C=O)]')
b.manage_atom_rules('Si', add=True)
b.process(remove_stereo=False)
b.save_results('test')

# 调用3：cas获取smiles后处理
from molppc import MolecularProcessor

c = MolecularProcessor()
c.load_data(input_data="Rabbit_embryo_loss.csv",
            cas_col='CAS Number', target_col='Rabbit_embryo_loss')
c.config_web_request(source='pubchem', max_workers=4)
c.add_web_request(cas=False, iupac=False, smiles=True)
c.process()
c.config_deduplicator(data_type='continuous')
c.deduplicate_data()
c.add_web_request(cas=False, iupac=True)
c.save_results('test2')

# 调用4：只去重smiles，用于预训练
from molppc import MolecularProcessor

d = MolecularProcessor()
d.load_data(input_data="ChEMBL_part3_.csv",
            smiles_col='Smiles')
d.process()
d.config_deduplicator()
d.deduplicate_data()
d.save_results('test3')

# 其他：自定义异常值检测逻辑
from molppc import MolecularProcessor
import pandas as pd
import numpy as np


def custom_outlier_filter(values: pd.Series):
    clean = values[(values >= -3.5) & (values <= 3.5)]
    if len(clean) > 1:
        return np.mean(clean), "custom_range"
    if len(clean) == 0:
        return np.mean(values), "all_custom_range"
    return clean, "custom_range"


e = MolecularProcessor()
e.load_data(input_data="logkOW_test.csv",
            smiles_col='Smiles', target_col='Value')
e.process()
e.config_deduplicator(custom_method=custom_outlier_filter, data_type='continuous')
e.deduplicate_data()
e.save_results('test4')
