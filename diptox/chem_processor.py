# diptox/chem_processor.py
from contextlib import redirect_stderr
from importlib import resources
import io
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, SaltRemover, rdmolops
from typing import List, Tuple, Optional, Callable
from .logger import log_manager
logger = log_manager.get_logger(__name__)


class ChemistryProcessor:
    """Handles all chemistry-related operations"""

    def __init__(self):
        self.remover = SaltRemover.SaltRemover()
        self._default_salts()
        self._solvents = self._default_solvents()
        self._neutralization_rules = self._default_neutralization_rules()
        self._valid = {'Br', 'C', 'Cl', 'F', 'H', 'I', 'N', 'O', 'P', 'S', 'Si', 'As', 'Se', 'Te', 'At'}
        self._custom_salts = []
        self._removed_salts = []
        self._custom_solvents = []
        self._removed_solvents = []

    @staticmethod
    def _default_neutralization_rules() -> List[Tuple[str, str]]:
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
            ('[$([N-;X2]C#N)]', 'N'),
            ('[c-H1]', '[CH2]'),
            ('[$([O-][N]C=O)]', 'O')
        ]

    def _default_salts(self):
        """Default salts."""
        salt_mols = []
        with resources.open_text("diptox", "salts.smi") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    name, smarts = line.split("\t")
                    mol = Chem.MolFromSmarts(smarts)
                    if mol:
                        salt_mols.append(mol)
        self.remover.salts = salt_mols

    @staticmethod
    def _default_solvents():
        """Default solvents."""
        solvent_mols = []
        with resources.open_text("diptox", "solvents.smi") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    name, smarts = line.split("\t")
                    mol = Chem.MolFromSmiles(smarts)
                    if mol:
                        solvent_mols.append(mol)
        return solvent_mols

    def _get_effective_salts(self) -> List[str]:
        """Get the current list of effective salts after additions and removals."""
        # Combine default and custom salts
        seen_smarts = set()
        combined_unique_mols = []

        for mol in self.remover.salts:
            smarts = Chem.MolToSmarts(mol)
            if smarts not in seen_smarts:
                seen_smarts.add(smarts)
                combined_unique_mols.append(mol)

        for mol in self._custom_salts:
            smarts = Chem.MolToSmarts(mol)
            if smarts not in seen_smarts:
                seen_smarts.add(smarts)
                combined_unique_mols.append(mol)

        removed_smarts = {Chem.MolToSmarts(m) for m in self._removed_salts}
        effective = [mol for mol in combined_unique_mols if Chem.MolToSmarts(mol) not in removed_smarts]

        return effective

    def _get_effective_solvents(self):
        """Get the current list of effective solvents after additions and removals."""
        seen_smiles = set()
        combined_unique_mols = []

        for mol in self._solvents:
            smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
            if smiles not in seen_smiles:
                seen_smiles.add(smiles)
                combined_unique_mols.append(mol)

        for mol in self._custom_solvents:
            smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
            if smiles not in seen_smiles:
                seen_smiles.add(smiles)
                combined_unique_mols.append(mol)

        removed_smiles = {Chem.MolToSmiles(m, isomericSmiles=True, canonical=True) for m in self._removed_solvents}
        effective = [mol for mol in combined_unique_mols if
                     Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True) not in removed_smiles]

        return effective

    def add_neutralization_rule(self, reactant: str, product: str):
        """
        Add a new neutralization rule to the list, ensuring the rule is valid and there are no conflicts.
        :param reactant: SMARTS string for the reactant.
        :param product: SMILES string for the product.
        """
        try:
            if not reactant or not product: return False
            patt = Chem.MolFromSmarts(reactant)
            if not patt:
                logger.warning(f"Invalid SMARTS pattern: {reactant}")
                return False
            repl = Chem.MolFromSmiles(product, sanitize=False)
            if not repl:
                logger.warning(f"Invalid SMILES string: {product}")
                return False
        except Exception as e:
            logger.error("Rule validation error", exc_info=True)
            return False

        for existing_reactant, existing_product in self._neutralization_rules:
            if existing_reactant == reactant:
                if existing_product != product:
                    logger.warning(
                        f"Rule conflict: SMARTS pattern '{reactant}' already exists, "
                        f"but the SMILES differs. Using the user-provided SMILES '{product}'")
                    self.remove_neutralization_rule(reactant)
                    self._neutralization_rules.append((reactant, product))
                    logger.info(f"Rule updated: {reactant} -> {product}")
                    return True
                else:
                    logger.info(f"Rule already exists, no need to add: {reactant} -> {product}")
                    return True

        self._neutralization_rules.append((reactant, product))
        logger.info(f"New rule added: {reactant} -> {product}")

        return True

    def remove_neutralization_rule(self, reactant: str):
        """Remove a matching rule from the neutralization rule list."""
        initial_len = len(self._neutralization_rules)
        self._neutralization_rules = [
            rule for rule in self._neutralization_rules if rule[0] != reactant
        ]
        if len(self._neutralization_rules) < initial_len:
            logger.info(f"Rule removed: {reactant}")
            return True
        else:
            logger.warning(f"No matching rule found: {reactant}")
            return False

    def add_effective_atom(self, atom):
        """Add a new atom symbol to the list of valid atoms."""
        if not isinstance(atom, str):
            return False
        try:
            if Chem.MolFromSmiles(f"[{atom}]") is None:
                logger.warning(f"Invalid atom symbol: {atom}")
                return False
        except:
            return False

        self._valid.add(atom)
        return True

    def delete_effective_atom(self, atom):
        """Remove an atom symbol from the list of valid atoms."""
        if atom in self._valid:
            self._valid.remove(atom)
            logger.info(f"Atom {atom} has been removed.")
            return True
        logger.info(f"Atom {atom} not found.")
        return False

    def add_default_salt(self, smarts: str):
        mol = Chem.MolFromSmarts(smarts)
        if mol:
            canon_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
            existing_smiles = {
                Chem.MolToSmiles(m, isomericSmiles=True, canonical=True)
                for m in self._custom_salts
            }
            if canon_smiles not in existing_smiles:
                self._custom_salts.append(mol)
                logger.info(f"Default salt added: {smarts}")
            return True
        else:
            logger.warning(f"Invalid SMILES: {smarts}")
            return False

    def remove_default_salt(self, smarts: str):
        target = Chem.MolFromSmarts(smarts)
        if not target:
            logger.warning(f"Invalid SMILES: {smarts}")
            return False
        target_smiles = Chem.MolToSmiles(target, isomericSmiles=True, canonical=True)
        existing_smiles = {
            Chem.MolToSmiles(m, isomericSmiles=True, canonical=True)
            for m in self._removed_salts
        }
        if target_smiles not in existing_smiles:
            self._removed_salts.append(target)
            logger.info(f"Default salt removed globally: {smarts}")
        return True

    def add_default_solvents(self, smarts: str):
        mol = Chem.MolFromSmarts(smarts)
        if mol:
            canon_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
            existing_smiles = {
                Chem.MolToSmiles(m, isomericSmiles=True, canonical=True)
                for m in self._custom_solvents
            }
            if canon_smiles not in existing_smiles:
                self._custom_solvents.append(mol)
                logger.info(f"Default solvent added: {smarts}")
                return True
        else:
            logger.warning(f"Invalid SMILES: {smarts}")
            return False

    def remove_default_solvents(self, smarts: str):
        target = Chem.MolFromSmarts(smarts)
        if not target:
            logger.warning(f"Invalid SMILES: {smarts}")
            return False
        target_smiles = Chem.MolToSmiles(target, isomericSmiles=True, canonical=True)
        existing_smiles = {
            Chem.MolToSmiles(m, isomericSmiles=True, canonical=True)
            for m in self._removed_solvents
        }
        if target_smiles not in existing_smiles:
            self._removed_solvents.append(target)
            logger.info(f"Default solvent removed globally: {smarts}")
        return True

    @staticmethod
    def mol_to_inchi(mol: Chem.Mol) -> Optional[str]:
        """Convert a molecule to InChI string using RDKit locally."""
        if mol is None:
            return None
        lg = RDLogger.logger()
        try:
            lg.setLevel(RDLogger.CRITICAL)
            inchi = Chem.MolToInchi(mol)
            return inchi
        except Exception as e:
            logger.warning(f"Failed to generate InChI: {str(e)}")
            return None
        finally:
            lg.setLevel(RDLogger.INFO)

    @staticmethod
    def CombineFragments(fragments: list[Chem.Mol]) -> Chem.Mol:
        """Combine multiple fragments into a single molecule."""
        smiles = [Chem.MolToSmiles(f, isomericSmiles=True, canonical=True) for f in fragments]
        return Chem.MolFromSmiles('.'.join(smiles))

    @staticmethod
    def smiles_to_mol(smiles: str, sanitize: bool = True) -> Optional[Chem.Mol]:
        """
        Convert a SMILES string to a molecular object.
        :param sanitize: Whether to perform chemical validation.
        """
        if smiles is None:
            return None
        with redirect_stderr(io.StringIO()):
            mol = Chem.MolFromSmiles(str(smiles), sanitize=sanitize)
        if mol is None:
            logger.warning(f"Invalid SMILES: {smiles}")
        return mol

    @staticmethod
    def standardize_smiles(mol: Chem.Mol, canonical: bool = True) -> Optional[str]:
        """
        Generate standardized SMILES.
        :param canonical: Whether to generate canonical form.
        """
        return Chem.MolToSmiles(mol, canonical=canonical)

    @staticmethod
    def remove_isotopes(mol: Chem.Mol) -> Chem.Mol:
        """
        Remove isotope information from all atoms in a molecule.
        This sets the isotope property of each atom to 0.
        """
        for atom in mol.GetAtoms():
            if atom.GetIsotope():
                atom.SetIsotope(0)
        return mol

    @staticmethod
    def reject_radicals(mol: Chem.Mol) -> Optional[Chem.Mol]:
        """If the molecule contains free radical electrons, it will be rejected"""
        if mol is None:
            return None
        for a in mol.GetAtoms():
            if a.GetNumRadicalElectrons() > 0:
                logger.warning(
                    f"Radical detected on atom {a.GetIdx()} "
                    f"({a.GetSymbol()}) in {Chem.MolToSmiles(mol)}. Molecule rejected."
                )
                return None
        return mol

    @staticmethod
    def remove_stereochemistry(mol: Chem.Mol) -> Chem.Mol:
        """Remove stereochemistry information"""
        Chem.RemoveStereochemistry(mol)
        return mol

    @staticmethod
    def remove_hydrogens(mol: Chem.Mol) -> Chem.Mol:
        """Remove hydrogen atoms"""
        return Chem.RemoveHs(mol)

    def neutralize_charges(self, mol: Chem.Mol, reject_non_neutral: bool = False) -> Optional[Chem.Mol]:
        """Charge neutralization processing and optionally reject non-neutral molecules."""
        for reactant, product in self._neutralization_rules:
            patt = Chem.MolFromSmarts(reactant)
            repl = Chem.MolFromSmiles(product, False)
            while mol.HasSubstructMatch(patt):
                mol = Chem.ReplaceSubstructs(mol, patt, repl)[0]
        Chem.SanitizeMol(mol)

        if reject_non_neutral:
            charge = rdmolops.GetFormalCharge(mol)
            if charge != 0:
                logger.info(f"Non-neutral molecule ({charge}) {Chem.MolToSmiles(mol)} rejected.")
                return None
        return mol

    def remove_salts(self, mol: Chem.Mol) -> Optional[Chem.Mol]:
        """Remove salts."""
        self.remover.salts = self._get_effective_salts()

        if mol.GetNumAtoms() >= 2 and len(Chem.GetMolFrags(mol)) == 1:
            return mol

        stripped = self.remover.StripMol(mol)
        if stripped.GetNumAtoms() == 0:
            frags = Chem.GetMolFrags(mol, asMols=True)
            best_frag = sorted(
                frags,
                key=lambda f: (
                    f.GetNumHeavyAtoms(),
                    sum(1 for a in f.GetAtoms() if a.GetAtomicNum() == 6),
                    Chem.MolToSmiles(f, canonical=True)
                ),
                reverse=True
            )[0]
            return best_frag

        fragments = Chem.GetMolFrags(stripped, asMols=True)
        if len(fragments) == 1:
            return stripped
        first_smiles = Chem.MolToSmiles(fragments[0], canonical=True)
        all_identical = all(
            Chem.MolToSmiles(frag, canonical=True) == first_smiles
            for frag in fragments[1:]
        )
        if all_identical:
            return fragments[0]
        else:
            combined = fragments[0]
            for frag in fragments[1:]:
                combined = self.CombineFragments([combined, frag])
            return combined

    def remove_solvents(self, mol: Chem.Mol) -> Optional[Chem.Mol]:
        """Remove solvents."""
        self._solvents = self._get_effective_solvents()

        fragments = Chem.GetMolFrags(mol, asMols=True)
        if len(fragments) == 1:
            return mol
        if len(fragments) == 0:
            return None

        non_solvent = [
            frag for frag in fragments
            if not any(frag.HasSubstructMatch(solvent) and solvent.HasSubstructMatch(frag) for solvent in self._solvents)
        ]
        if not non_solvent:
            best_frag = sorted(
                fragments,
                key=lambda f: (
                    f.GetNumHeavyAtoms(),
                    sum(1 for a in f.GetAtoms() if a.GetAtomicNum() == 6),
                    Chem.MolToSmiles(f, canonical=True)
                ),
                reverse=True
            )[0]
            return best_frag

        first_smiles = Chem.MolToSmiles(non_solvent[0], canonical=True)
        all_identical = all(
            Chem.MolToSmiles(frag, canonical=True) == first_smiles
            for frag in non_solvent[1:]
        )
        if all_identical:
            return non_solvent[0]
        else:
            combined = non_solvent[0]
            for frag in non_solvent[1:]:
                combined = self.CombineFragments([combined, frag])
            return combined

    @staticmethod
    def remove_mixtures(mol: Chem.Mol,
                        hac_threshold: int = 3,
                        keep_largest: bool = True,) -> Optional[Chem.Mol]:
        """Remove mixtures.
        :param keep_largest: Whether to keep the largest molecules.
        :param hac_threshold: Threshold of hac molecules."""
        fragments = list(rdmolops.GetMolFrags(mol, asMols=True))
        if len(fragments) > 1:
            fragments = [f for f in fragments if f.GetNumHeavyAtoms() > hac_threshold]
        if len(fragments) > 1:
            logger.warning(
                f"{Chem.MolToSmiles(mol)} contains >1 fragment with >" + str(hac_threshold) + " heavy atoms")
            return max(fragments, key=lambda x: x.GetNumHeavyAtoms()) if keep_largest else None
        elif len(fragments) == 0:
            logger.warning(
                f"{Chem.MolToSmiles(mol)} contains no fragments with >" + str(hac_threshold) + " heavy atoms")
            return None
        else:
            return fragments[0]

    @staticmethod
    def remove_inorganic(mol: Chem.Mol) -> Optional[Chem.Mol]:
        """Remove inorganic atoms"""
        has_carbon = any(atom.GetSymbol() == 'C' for atom in mol.GetAtoms())
        if not has_carbon:
            return None

        inorganic_smarts = [
            '[#6]#[#7]', '[#8]=[#6]=[#8]', '[#6]#[#8]', '[#8]=[#6]=[#16]', '[#8]=[#6]=[#6]=[#6]=[#8]',
            '[#8]~[#6](=[#8])~[#8]', '[#16]=[#6]=[#16]', '[#7]=[#6]=[#8]', '[#8]-[#6]#[#7]',
            '[#16]-[#6]#[#7]', '[#7]=[#6]=[#16]', '[#34]-[#6]#[#7]', '[F,Cl,Br,I]-[#6]#[#7]',
            '[#7]#[#6]-[#6]#[#7]', '[Cl]-[#6](=[#8])-[Cl]', '[Cl]-[#6](=[#16])-[Cl]', '[#6]',
        ]
        for smarts in inorganic_smarts:
            pattern = Chem.MolFromSmarts(smarts)
            if not pattern:
                continue
            if mol.GetNumAtoms() == pattern.GetNumAtoms() and mol.HasSubstructMatch(pattern):
                return None

        return mol

    def effective_atom(self, mol: Chem.Mol, strict: bool = False) -> Optional[Chem.Mol]:
        """
        Check if all atoms in the molecule are valid according to the list of valid atoms.
        :param mol: The molecular object to check.
        :param strict: If True, reject the entire molecule if any invalid atoms are found.
                       If False, attempt to remove only terminal invalid atoms.
        """
        invalid_elements = {atom.GetSymbol() for atom in mol.GetAtoms()} - self._valid
        if not invalid_elements:
            return mol
        original_smiles = Chem.MolToSmiles(mol)
        if strict:
            logger.error(
                f"Strict check: Invalid atoms {invalid_elements} detected in molecule {original_smiles}. Molecule removed."
            )
            return None

        rw_mol = Chem.RWMol(mol)
        atoms_to_remove_indices = []

        for atom in rw_mol.GetAtoms():
            if atom.GetSymbol() in invalid_elements:
                if atom.GetDegree() > 1:
                    logger.warning(
                        f"Cannot remove non-terminal invalid atom '{atom.GetSymbol()}' from {original_smiles}. Molecule rejected."
                    )
                    return None
                atoms_to_remove_indices.append(atom.GetIdx())

        for idx in sorted(atoms_to_remove_indices, reverse=True):
            rw_mol.RemoveAtom(idx)

        if rw_mol.GetNumAtoms() == 0:
            logger.warning(f"Molecule {original_smiles} became empty after removing invalid atoms. Molecule rejected.")
            return None

        try:
            cleaned_mol = rw_mol.GetMol()
            Chem.SanitizeMol(cleaned_mol)
            logger.info(
                f"Removed invalid atoms {invalid_elements} from {original_smiles}, "
                f"resulting in {Chem.MolToSmiles(cleaned_mol)}."
            )
            return cleaned_mol
        except Exception as e:
            logger.error(
                f"Sanitization failed for {original_smiles} after removing atoms: {e}. Molecule rejected."
            )
            return None

    def display_current_rules(self) -> None:
        """Prints a summary of the currently active chemical processing rules."""
        print("--- Current Chemical Processing Rules ---")

        # Valid Atoms
        print("\n[+] Valid Atoms:")
        print(f"    {', '.join(sorted(list(self._valid)))}")

        # Neutralization Rules
        print("\n[+] Neutralization Rules (Reactant SMARTS -> Product SMILES):")
        for reactant, product in self._neutralization_rules:
            print(f"    - {reactant} -> {product}")

        # Salts
        salts = self._get_effective_salts()
        print(f"\n[+] Effective Salts ({len(salts)} total):")
        for salt_mol in salts:
            print(f"    - {Chem.MolToSmarts(salt_mol)}")

        # Solvents
        solvents = self._get_effective_solvents()
        print(f"\n[+] Effective Solvents ({len(solvents)} total):")
        for solvent_mol in solvents:
            print(f"    - {Chem.MolToSmiles(solvent_mol)}")

        print("\n--- End of Rules ---")

    def get_current_rules_dict(self) -> dict:
        """Return current rules as a dictionary for GUI display."""
        current_salts = []
        for mol in self._get_effective_salts():
            try:
                current_salts.append(Chem.MolToSmarts(mol))
            except:
                current_salts.append("Invalid Salt Pattern")

        current_solvents = []
        for mol in self._get_effective_solvents():
            try:
                current_solvents.append(Chem.MolToSmiles(mol, isomericSmiles=True))
            except:
                current_solvents.append("Invalid Solvent Pattern")

        return {
            "atoms": sorted(list(self._valid)),
            "salts": current_salts,
            "solvents": current_solvents,
            "neutralization": self._neutralization_rules
        }

    @staticmethod
    def validate_atom_count(mol: Optional[Chem.Mol],
                            min_heavy_atoms: Optional[int] = None,
                            max_heavy_atoms: Optional[int] = None,
                            min_total_atoms: Optional[int] = None,
                            max_total_atoms: Optional[int] = None) -> bool:
        """
        Check if a molecule is within the specified atom count limits.
        :param mol: The RDKit molecule object to check.
        :param min_heavy_atoms: Minimum number of heavy atoms (inclusive).
        :param max_heavy_atoms: Maximum number of heavy atoms (inclusive).
        :param min_total_atoms: Minimum number of total atoms (inclusive).
        :param max_total_atoms: Maximum number of total atoms (inclusive).
        :return: True if the molecule meets all criteria, False otherwise.
        """
        if mol is None:
            return False

        if min_heavy_atoms is not None and mol.GetNumHeavyAtoms() < min_heavy_atoms:
            return False

        if max_heavy_atoms is not None and mol.GetNumHeavyAtoms() > max_heavy_atoms:
            return False

        if min_total_atoms is not None or max_total_atoms is not None:
            try:
                mol.UpdatePropertyCache(strict=False)

                mol_with_hs = Chem.AddHs(mol)
                total_atoms = mol_with_hs.GetNumAtoms()

                if min_total_atoms is not None and total_atoms < min_total_atoms:
                    return False

                if max_total_atoms is not None and total_atoms > max_total_atoms:
                    return False
            except Exception:
                return False

        return True

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
            lambda mol: processor_instance.remove_salts(mol),
            processor_instance.remove_inorganic,
            processor_instance.neutralize_charges,
            processor_instance.effective_atom,
        )
