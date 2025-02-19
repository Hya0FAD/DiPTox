# molppc/chem_processor.py
from contextlib import redirect_stderr
import io
from rdkit import Chem
from rdkit.Chem import AllChem, SaltRemover, rdmolops
from typing import List, Tuple, Optional, Callable
from .logger import log_manager
logger = log_manager.get_logger(__name__)


class ChemistryProcessor:
    """Handles all chemistry-related operations"""

    def __init__(self):
        self.remover = SaltRemover.SaltRemover()
        self._neutralization_rules = self._default_rules()
        self._valid = {'Br', 'C', 'Cl', 'F', 'H', 'I', 'N', 'O', 'P', 'S'}

    @staticmethod
    def _default_rules() -> List[Tuple[str, str]]:
        """Default charge neutralization rules"""
        return [
            ('[n+;H]', 'n'),
            ('[N+;!H0]', 'N'),
            ('[$([O-]);!$([O-][#7])]', 'O'),
            ('[S-;X1]', 'S'),
            ('[$([N-;X2]S(=O)=O)]', 'N'),
            ('[$([N-;X2][C,N]=C)]', 'N'),
            ('[n-]', '[nH]'),
            ('[$([S-]=O)]', 'S'),
            ('[$([N-]C=O)]', 'N'),
        ]

    def add_neutralization_rule(self, reactant: str, product: str):
        """
        Add a new neutralization rule to the list, ensuring the rule is valid and there are no conflicts.
        :param reactant: SMARTS string for the reactant.
        :param product: SMILES string for the product.
        """
        try:
            patt = Chem.MolFromSmarts(reactant)
            if not patt:
                logger.warning(f"Invalid SMARTS pattern: {reactant}")
                return

            repl = Chem.MolFromSmiles(product, sanitize=False)
            if not repl:
                logger.warning(f"Invalid SMILES string: {product}")
                return

        except Exception as e:
            logger.error("Rule validation error", exc_info=True)
            return

        for existing_reactant, existing_product in self._neutralization_rules:
            if existing_reactant == reactant:
                if existing_product != product:
                    logger.warning(
                        f"Rule conflict: SMARTS pattern '{reactant}' already exists, "
                        f"but the SMILES differs. Using the user-provided SMILES '{product}'")
                    self.remove_neutralization_rule(reactant)
                    self._neutralization_rules.append((reactant, product))
                    logger.info(f"Rule updated: {reactant} -> {product}")
                    return
                else:
                    logger.info(f"Rule already exists, no need to add: {reactant} -> {product}")
                    return

        self._neutralization_rules.append((reactant, product))
        logger.info(f"New rule added: {reactant} -> {product}")

    def remove_neutralization_rule(self, reactant: str):
        """
        Remove a matching rule from the neutralization rule list.
        """
        initial_len = len(self._neutralization_rules)
        self._neutralization_rules = [
            rule for rule in self._neutralization_rules if rule[0] != reactant
        ]
        if len(self._neutralization_rules) < initial_len:
            logger.info(f"Rule removed: {reactant}")
        else:
            logger.warning(f"No matching rule found: {reactant}")

    def add_effective_atom(self, atom):
        """
        Add a new atom symbol to the list of valid atoms.
        """
        if isinstance(atom, str):
            self._valid.add(atom)
        else:
            logger.info(f"Invalid atom format: {atom}. It should be a string representing an atom symbol.")

    def delete_effective_atom(self, atom):
        """
        Remove an atom symbol from the list of valid atoms.
        """
        if isinstance(atom, str):
            if atom in self._valid:
                self._valid.remove(atom)
                logger.info(f"Atom {atom} has been removed from valid atoms.")
            else:
                logger.info(f"Atom {atom} is not in the valid atoms list. No changes made.")
        else:
            logger.info(f"Invalid atom format: {atom}. It should be a string representing an atom symbol.")

    @staticmethod
    def smiles_to_mol(smiles: str, sanitize: bool = True) -> Optional[Chem.Mol]:
        """
        Convert a SMILES string to a molecular object.
        :param sanitize: Whether to perform chemical validation.
        """
        with redirect_stderr(io.StringIO()):
            mol = Chem.MolFromSmiles(str(smiles), sanitize=sanitize)
        if mol is None:
            logger.warning(f"Invalid SMILES: {smiles}")
        return mol

    @staticmethod
    def standardize_smiles(mol: Chem.Mol,
                           isomeric: bool = False,
                           canonical: bool = True) -> Optional[str]:
        """
        Generate standardized SMILES.
        :param isomeric: Whether to retain stereochemistry.
        :param canonical: Whether to generate canonical form.
        """
        return Chem.MolToSmiles(mol, isomericSmiles=isomeric, canonical=canonical)

    @staticmethod
    def remove_stereochemistry(mol: Chem.Mol) -> Chem.Mol:
        """Remove stereochemistry information"""
        Chem.RemoveStereochemistry(mol)
        return mol

    @staticmethod
    def remove_hydrogens(mol: Chem.Mol) -> Chem.Mol:
        """Remove hydrogen atoms"""
        return Chem.RemoveHs(mol)

    def neutralize_charges(self, mol: Chem.Mol) -> Chem.Mol:
        """Charge neutralization processing"""
        for reactant, product in self._neutralization_rules:
            patt = Chem.MolFromSmarts(reactant)
            repl = Chem.MolFromSmiles(product, False)
            while mol.HasSubstructMatch(patt):
                mol = Chem.ReplaceSubstructs(mol, patt, repl)[0]
        return mol

    def remove_salts(self, mol: Chem.Mol,
                     hac_threshold: int = 3,
                     keep_largest: bool = False) -> Optional[Chem.Mol]:
        """
        Remove salts/solvents.
        :param hac_threshold: Retain fragments with more than this number of heavy atoms.
        :param keep_largest: Whether to keep the largest fragment instead of the first valid fragment.
        """
        stripped = self.remover.StripMol(mol)
        fragments = list(rdmolops.GetMolFrags(stripped, asMols=True))
        if len(fragments) > 1:
            fragments = [f for f in fragments if f.GetNumAtoms() > hac_threshold]
        if len(fragments) > 1:
            logger.warning(f"{Chem.MolToSmiles(mol)} contains >1 fragment with >" + str(hac_threshold) + " heavy atoms")
            return None
        elif len(fragments) == 0:
            logger.warning(f"{Chem.MolToSmiles(mol)} contains no fragments with >" + str(hac_threshold) + " heavy atoms")
            return None
        else:
            return max(fragments, key=lambda x: x.GetNumAtoms()) if keep_largest else fragments[0]

    def effective_atom(self, mol: Chem.Mol) -> Optional[Chem.Mol]:
        """
        Check if all atoms in the molecule are valid according to the list of valid atoms.
        :param mol: The molecular object to check.
        """
        try:
            elements = [atom.GetSymbol() for atom in mol.GetAtoms()]
            if all(atom in self._valid for atom in elements):
                return mol
            else:
                invalid_atoms = set(elements) - self._valid
                logger.error(
                    f"Invalid atoms detected: {invalid_atoms} "
                    f"in molecule {Chem.MolToSmiles(mol)}"
                )
        except Exception as e:
            logger.exception("Atom validation failed")
            return None

    @classmethod
    def create_pipeline(cls, *processors: Callable[[Chem.Mol], Chem.Mol]):
        """Create a custom processing pipeline"""

        def pipeline(mol: Chem.Mol):
            for processor in processors:
                if mol is None:  # Allow the process to terminate early if any step returns None
                    break
                mol = processor(mol)
            return mol

        return pipeline

    @classmethod
    def default_standardization(cls, processor_instance) -> Callable:
        """Get the default standardization pipeline"""
        return cls.create_pipeline(
            processor_instance.remove_hydrogens,
            processor_instance.remove_stereochemistry,
            processor_instance.neutralize_charges,
            lambda mol: processor_instance.remove_salts(mol),
            processor_instance.effective_atom,
        )
