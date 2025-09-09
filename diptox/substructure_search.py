# diptox/substructure_search.py
from rdkit import Chem
import pandas as pd
from typing import Tuple, List, Dict
from .logger import log_manager
logger = log_manager.get_logger(__name__)


class SubstructureSearcher:

    def __init__(self,
                 df: pd.DataFrame,
                 smiles_col: str,
                 processed_mols: List[Chem.Mol] = None):
        self.df = df
        self.smiles_col = smiles_col
        self.processed_mols = processed_mols
        self._raw_mols = self._load_raw_molecules() if processed_mols is None else None

    def _load_raw_molecules(self) -> List[Tuple[int, Chem.Mol]]:
        """Load the unprocessed molecule and keep the original index"""
        return [
            (idx, Chem.MolFromSmiles(smi))
            for idx, smi in self.df[self.smiles_col].items()
            if pd.notna(smi)
        ]

    @staticmethod
    def _parse_query(pattern: str, is_smarts: bool) -> Chem.Mol:
        """Parsing query mode (SMILES/SMARTS support)"""
        mol = Chem.MolFromSmarts(pattern) if is_smarts else Chem.MolFromSmiles(pattern)
        if not mol:
            raise ValueError(f"Invalid {'SMARTS' if is_smarts else 'SMILES'}: {pattern}")
        return mol

    def search(self,
               query_pattern: str,
               is_smarts: bool = False,
               return_mols: bool = False) -> Dict:
        """
        Perform a substructure search
        :return: {
            'count': Number of matches,
            'matches': [(raw data index, SMILES or Mol object),...] ,
            'query_type': 'SMARTS'/'SMILES'
        }
        """
        try:
            query_mol = self._parse_query(query_pattern, is_smarts)
        except Exception as e:
            logger.error(str(e))
            return {'count': 0, 'matches': [], 'query_type': 'Invalid'}

        if self.processed_mols:
            source = zip(self.df.index, self.processed_mols)
        else:
            source = [(idx, mol) for idx, mol in self._raw_mols if mol]

        matches = []
        for idx, mol in source:
            if mol.HasSubstructMatch(query_mol, useChirality=True):
                result = (idx, mol if return_mols else Chem.MolToSmiles(mol, isomericSmiles=True))
                matches.append(result)
        logger.info(f"{len(matches)} matches were found with {query_pattern}")
        return {
            'count': len(matches),
            'matches': matches,
            'query_type': 'SMARTS' if is_smarts else 'SMILES'
        }
