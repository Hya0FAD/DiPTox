# diptox/chem_processor.py
from contextlib import redirect_stderr
from importlib import resources
import io
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, SaltRemover, rdmolops
from typing import List, Tuple, Optional, Callable
from .element_policy import REMOVABLE_ALKALI_METALS, is_metal
from .logger import log_manager
logger = log_manager.get_logger(__name__)


class AmbiguousFragmentError(ValueError):
    """The largest-fragment policy cannot choose a unique molecular identity."""

    def __init__(self, heavy_atom_count: int, candidates: List[str]):
        self.heavy_atom_count = heavy_atom_count
        self.candidates = tuple(sorted(candidates))
        super().__init__(
            f"Ambiguous largest fragment: {len(self.candidates)} distinct components "
            f"have {heavy_atom_count} heavy atoms ({'; '.join(self.candidates)})"
        )


class ChemistryProcessor:
    """Handles all chemistry-related operations"""

    def __init__(self):
        self.remover = SaltRemover.SaltRemover()
        self._default_salts()
        self._solvents = self._default_solvents()
        self._neutralization_rules = self._default_neutralization_rules()
        self._valid = {'B', 'Br', 'C', 'Cl', 'F', 'H', 'I', 'N', 'O', 'P', 'S', 'Si', 'As', 'Se', 'Te', 'At'}
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
            ('[$([n-]1nnnc1),$([n-]1nncn1)]', '[nH]'),
            ('[$([N-]C=O)]', 'N'),
            ('[$([N-;X2]C#N)]', 'N'),
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
    def CombineFragments(fragments: List[Chem.Mol]) -> Chem.Mol:
        """Combine multiple fragments into a single molecule."""
        combined = Chem.Mol()
        for fragment in fragments:
            combined = Chem.CombineMols(combined, fragment)
        return combined

    @staticmethod
    def smiles_to_mol(smiles: str, sanitize: bool = True,
                      remove_hs: bool = True) -> Optional[Chem.Mol]:
        """
        Convert a SMILES string to a molecular object.
        :param sanitize: Whether to perform chemical validation.
        :param remove_hs: Whether the parser should remove ordinary explicit hydrogen atoms.
        """
        if not isinstance(smiles, str):
            return None
        parameters = Chem.SmilesParserParams()
        parameters.sanitize = sanitize
        parameters.removeHs = remove_hs
        with redirect_stderr(io.StringIO()):
            mol = Chem.MolFromSmiles(smiles, parameters)
        if mol is None:
            logger.warning(f"Invalid SMILES: {smiles}")
        return mol

    @staticmethod
    def standardize_smiles(mol: Chem.Mol, canonical: bool = True) -> Optional[str]:
        """
        Generate standardized SMILES without atom-map annotations, leaving the input intact.
        :param canonical: Whether to generate canonical form.
        """
        unannotated = Chem.Mol(mol)
        for atom in unannotated.GetAtoms():
            atom.SetAtomMapNum(0)
        return Chem.MolToSmiles(unannotated, canonical=canonical)

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
        """Reject radicals except the known nitric oxide and aminoxyl motifs."""
        if mol is None:
            return None
        for a in mol.GetAtoms():
            radical_electrons = a.GetNumRadicalElectrons()
            if radical_electrons == 0:
                continue

            neighbors = a.GetNeighbors()
            allowed = False
            if radical_electrons == 1 and len(neighbors) == 1:
                neighbor = neighbors[0]
                bond = mol.GetBondBetweenAtoms(a.GetIdx(), neighbor.GetIdx())
                allowed = (
                    a.GetAtomicNum() == 7
                    and neighbor.GetAtomicNum() == 8
                    and bond.GetBondType() == Chem.BondType.DOUBLE
                ) or (
                    a.GetAtomicNum() == 8
                    and neighbor.GetAtomicNum() == 7
                    and bond.GetBondType() == Chem.BondType.SINGLE
                )

            if not allowed:
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

    @staticmethod
    def add_hydrogens(mol: Chem.Mol) -> Chem.Mol:
        """Expand implicit hydrogen counts to explicit atoms on a molecule copy."""
        return Chem.AddHs(mol)

    @staticmethod
    def reject_dummy_atoms(mol: Chem.Mol) -> Optional[Chem.Mol]:
        """Reject wildcard/dummy atoms because they do not define a complete molecule."""
        return None if any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms()) else mol

    @staticmethod
    def _replace_neutralization_site(mol: Chem.Mol, pattern: Chem.Mol,
                                    replacement: Chem.Mol,
                                    proton_transfer: bool) -> Chem.Mol:
        """Preserve atom identity for single-atom rules; retain general custom replacements."""
        # A fixed charge can leave only some acid/base sites neutralized. Choose
        # those sites in canonical graph order, independently of input atom order
        # and atom-map labels, without renumbering the retained source molecule.
        ranking_copy = Chem.Mol(mol)
        for atom in ranking_copy.GetAtoms():
            atom.SetAtomMapNum(0)
        ranks = Chem.CanonicalRankAtoms(ranking_copy, includeChirality=True,
                                       includeIsotopes=True)
        order = sorted(range(mol.GetNumAtoms()), key=lambda idx: ranks[idx])
        ordered = Chem.RenumberAtoms(mol, order)
        match = tuple(order[idx] for idx in ordered.GetSubstructMatch(pattern))
        if pattern.GetNumAtoms() != 1 or replacement.GetNumAtoms() != 1:
            return Chem.ReplaceSubstructs(ordered, pattern, replacement, replaceAll=False)[0]

        original = mol.GetAtomWithIdx(match[0])
        product = replacement.GetAtomWithIdx(0)
        if original.GetAtomicNum() != product.GetAtomicNum():
            return Chem.ReplaceSubstructs(ordered, pattern, replacement, replaceAll=False)[0]

        candidate = Chem.RWMol(mol)
        atom = candidate.GetAtomWithIdx(match[0])
        hydrogen_indices = []
        if proton_transfer:
            hydrogen_neighbors = [neighbor for neighbor in original.GetNeighbors()
                                  if neighbor.GetAtomicNum() == 1]
            target_hydrogens = (original.GetTotalNumHs(includeNeighbors=True)
                                + product.GetFormalCharge() - original.GetFormalCharge())
            remove_count = max(0, len(hydrogen_neighbors) - target_hydrogens)
            # Prefer unlabelled protons when the source explicitly represents every H.
            hydrogen_neighbors.sort(key=lambda h: (h.GetIsotope() != 0,
                                                   h.GetAtomMapNum() != 0,
                                                   h.GetIsotope(), h.GetAtomMapNum(),
                                                   ranks[h.GetIdx()]))
            hydrogen_indices = [h.GetIdx() for h in hydrogen_neighbors[:remove_count]]
            atom.SetNumExplicitHs(max(0, target_hydrogens
                                      - len(hydrogen_neighbors) + len(hydrogen_indices)))
            atom.SetNoImplicit(True)
        else:
            atom.SetNumExplicitHs(product.GetNumExplicitHs())
            atom.SetNoImplicit(product.GetNoImplicit())

        atom.SetFormalCharge(product.GetFormalCharge())
        atom.SetIsAromatic(product.GetIsAromatic())
        atom.SetNumRadicalElectrons(product.GetNumRadicalElectrons())
        if product.GetIsotope():
            atom.SetIsotope(product.GetIsotope())
        if product.GetAtomMapNum():
            atom.SetAtomMapNum(product.GetAtomMapNum())

        for hydrogen_idx in sorted(hydrogen_indices, reverse=True):
            neighbors = [neighbor.GetIdx() for neighbor in atom.GetNeighbors()]
            if (atom.GetChiralTag() in (Chem.ChiralType.CHI_TETRAHEDRAL_CW,
                                       Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
                    and (len(neighbors) - 1 - neighbors.index(hydrogen_idx)) % 2):
                atom.InvertChirality()
            candidate.RemoveAtom(hydrogen_idx)
        return candidate.GetMol()

    def neutralize_charges(self, mol: Chem.Mol, reject_non_neutral: bool = False) -> Optional[Chem.Mol]:
        """Neutralize proton-transfer sites while preserving fixed-charge compensation."""
        replacement_limit = max(100, 10 * mol.GetNumAtoms())
        replacements = 0
        compiled_rules = []
        neutralizable_atoms = set()
        default_rules = set(self._default_neutralization_rules())
        for reactant, product in self._neutralization_rules:
            patt = Chem.MolFromSmarts(reactant)
            repl = Chem.MolFromSmiles(product, False)
            compiled_rules.append((patt, repl, (reactant, product) in default_rules))
            for match in mol.GetSubstructMatches(patt):
                if match and mol.GetAtomWithIdx(match[0]).GetFormalCharge() != 0:
                    neutralizable_atoms.add(match[0])

        # Balanced fixed charges, such as a nitro group's N+/O-, do not need
        # compensation from proton-transfer sites elsewhere in the molecule.
        permanent_net_charge = sum(
            atom.GetFormalCharge() for atom in mol.GetAtoms()
            if atom.GetIdx() not in neutralizable_atoms
        )

        for patt, repl, proton_transfer in compiled_rules:
            seen_structures = {Chem.MolToSmiles(mol, canonical=True)}
            while mol.HasSubstructMatch(patt):
                if replacements >= replacement_limit:
                    logger.warning("Neutralization stopped: custom rules exceeded the replacement limit.")
                    return None
                old_charge = rdmolops.GetFormalCharge(mol)
                candidate = self._replace_neutralization_site(mol, patt, repl, proton_transfer)
                Chem.SanitizeMol(candidate)
                new_charge = rdmolops.GetFormalCharge(candidate)
                if permanent_net_charge and abs(new_charge) >= abs(old_charge):
                    break
                candidate_smiles = Chem.MolToSmiles(candidate, canonical=True)
                if candidate_smiles in seen_structures:
                    logger.warning("Neutralization stopped: a custom rule repeated a structure without progress.")
                    return None
                seen_structures.add(candidate_smiles)
                replacements += 1
                mol = candidate
        Chem.SanitizeMol(mol)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

        if reject_non_neutral and self.reject_non_neutral(mol) is None:
            return None
        return mol

    @staticmethod
    def reject_non_neutral(mol: Chem.Mol) -> Optional[Chem.Mol]:
        """Reject a molecule whose total formal charge is not zero."""
        charge = rdmolops.GetFormalCharge(mol)
        if charge != 0:
            logger.info(f"Non-neutral molecule ({charge}) {Chem.MolToSmiles(mol)} rejected.")
            return None
        return mol

    @staticmethod
    def _fragment_matches(fragment: Chem.Mol, patterns: List[Chem.Mol]) -> bool:
        """Match whole fragments, allowing only a narrow sulfoxide representation fallback."""
        # Match implicit-H dictionaries without changing explicitly retained source atoms.
        candidates = [fragment]
        if any(atom.GetAtomicNum() == 1 and atom.GetDegree() > 0
               and atom.GetIsotope() == 0 for atom in fragment.GetAtoms()):
            candidates.append(Chem.RemoveHs(fragment))

        def matches(candidate, pattern):
            return (candidate.GetNumAtoms() == pattern.GetNumAtoms()
                    and candidate.GetNumBonds() == pattern.GetNumBonds()
                    and candidate.HasSubstructMatch(pattern, useChirality=True))

        if any(matches(candidate, pattern) for pattern in patterns for candidate in candidates):
            return True
        return any(
            matches(view, pattern)
            for candidate in candidates
            for view in ChemistryProcessor._sulfoxide_matching_views(candidate)
            for pattern in patterns
        )

    @staticmethod
    def _sulfoxide_matching_views(fragment: Chem.Mol) -> List[Chem.Mol]:
        """Copy R-S(=O)-R / R-S+(-O-)-R views; never normalize SMARTS or retained molecules."""
        if fragment.HasQuery() or any(atom.GetNumRadicalElectrons() for atom in fragment.GetAtoms()):
            return []
        sites = []
        for sulfur in fragment.GetAtoms():
            if (sulfur.GetAtomicNum() != 16 or sulfur.GetIsAromatic()
                    or sulfur.GetDegree() != 3 or sulfur.GetTotalNumHs(includeNeighbors=True)):
                continue
            neighbors = list(sulfur.GetNeighbors())
            oxygens = [atom for atom in neighbors if atom.GetAtomicNum() == 8]
            carbons = [atom for atom in neighbors if atom.GetAtomicNum() == 6]
            if len(oxygens) != 1 or len(carbons) != 2:
                continue
            oxygen = oxygens[0]
            if oxygen.GetDegree() != 1 or oxygen.GetTotalNumHs(includeNeighbors=True):
                continue
            if any(fragment.GetBondBetweenAtoms(sulfur.GetIdx(), carbon.GetIdx()).GetBondType()
                   != Chem.BondType.SINGLE for carbon in carbons):
                continue
            bond = fragment.GetBondBetweenAtoms(sulfur.GetIdx(), oxygen.GetIdx())
            neutral = (sulfur.GetFormalCharge() == 0 and oxygen.GetFormalCharge() == 0
                       and bond.GetBondType() == Chem.BondType.DOUBLE)
            separated = (sulfur.GetFormalCharge() == 1 and oxygen.GetFormalCharge() == -1
                         and bond.GetBondType() == Chem.BondType.SINGLE)
            if neutral or separated:
                sites.append((sulfur.GetIdx(), oxygen.GetIdx(), separated))

        views = []
        for separated in (False, True):
            if not any(current != separated for _, _, current in sites):
                continue
            view = Chem.RWMol(fragment)
            for sulfur_idx, oxygen_idx, _ in sites:
                view.GetAtomWithIdx(sulfur_idx).SetFormalCharge(1 if separated else 0)
                view.GetAtomWithIdx(oxygen_idx).SetFormalCharge(-1 if separated else 0)
                view.GetBondBetweenAtoms(sulfur_idx, oxygen_idx).SetBondType(
                    Chem.BondType.SINGLE if separated else Chem.BondType.DOUBLE
                )
            molecule = view.GetMol()
            if Chem.SanitizeMol(molecule, catchErrors=True) == Chem.SanitizeFlags.SANITIZE_NONE:
                Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
                views.append(molecule)
        return views

    @staticmethod
    def _collapse_identical_fragments(fragments: List[Chem.Mol]) -> List[Chem.Mol]:
        """Retain one copy only when every fragment has the same map-free identity."""
        if len(fragments) <= 1:
            return fragments
        identity = ChemistryProcessor.standardize_smiles(fragments[0])
        if all(ChemistryProcessor.standardize_smiles(fragment) == identity
               for fragment in fragments[1:]):
            return fragments[:1]
        return fragments

    @staticmethod
    def collapse_identical_components(mol: Chem.Mol) -> Chem.Mol:
        """Fold an all-identical final mixture without partially deduplicating A.A.B."""
        fragments = list(Chem.GetMolFrags(mol, asMols=True))
        retained = ChemistryProcessor._collapse_identical_fragments(fragments)
        return retained[0] if len(fragments) > 1 and len(retained) == 1 else mol

    @staticmethod
    def _has_carbon(fragment: Chem.Mol) -> bool:
        return any(atom.GetAtomicNum() == 6 for atom in fragment.GetAtoms())

    @staticmethod
    def _contains_non_counterion_metal(fragment: Chem.Mol) -> bool:
        return any(
            is_metal(atom.GetAtomicNum()) and atom.GetAtomicNum() not in REMOVABLE_ALKALI_METALS
            for atom in fragment.GetAtoms()
        )

    @staticmethod
    def _contains_metal(fragment: Chem.Mol) -> bool:
        return any(is_metal(atom.GetAtomicNum()) for atom in fragment.GetAtoms())

    def remove_salts(self, mol: Chem.Mol) -> Optional[Chem.Mol]:
        """Remove counterions while protecting a possible organic parent fragment."""
        fragments = list(Chem.GetMolFrags(mol, asMols=True))
        if len(fragments) <= 1:
            return mol

        salts = self._get_effective_salts()
        solvents = self._get_effective_solvents()
        salt_flags = [self._fragment_matches(fragment, salts) for fragment in fragments]
        non_salts = [fragment for fragment, is_salt in zip(fragments, salt_flags) if not is_salt]
        carbon_salts = [
            fragment for fragment, is_salt in zip(fragments, salt_flags)
            if is_salt and self._has_carbon(fragment)
        ]
        protected_metals = [
            fragment for fragment, is_salt in zip(fragments, salt_flags)
            if is_salt and not self._has_carbon(fragment)
            and self._contains_non_counterion_metal(fragment)
        ]
        non_salt_parents = [
            fragment for fragment in non_salts
            if self._has_carbon(fragment) and not self._fragment_matches(fragment, solvents)
        ]

        if non_salt_parents:
            retained = non_salts + protected_metals
        elif carbon_salts:
            retained = carbon_salts + non_salts + protected_metals
        else:
            # Do not turn a salt/solvent mixture into a solvent-only success.
            retained = fragments

        return self.CombineFragments(self._collapse_identical_fragments(retained))

    def remove_solvents(self, mol: Chem.Mol) -> Optional[Chem.Mol]:
        """Remove solvents."""
        self._solvents = self._get_effective_solvents()

        fragments = Chem.GetMolFrags(mol, asMols=True)
        if len(fragments) == 1:
            return mol
        if len(fragments) == 0:
            return None

        non_solvent = [frag for frag in fragments if not self._fragment_matches(frag, self._solvents)]
        retained = self._collapse_identical_fragments(non_solvent or list(fragments))
        if len(retained) == 1:
            return retained[0]
        return self.CombineFragments(retained)

    @staticmethod
    def remove_mixtures(mol: Chem.Mol,
                        hac_threshold: int = 3,
                        keep_largest: bool = False,
                        mode: Optional[str] = None) -> Optional[Chem.Mol]:
        """Handle disconnected components by keeping, rejecting, or selecting the largest."""
        fragments = ChemistryProcessor._collapse_identical_fragments(
            list(rdmolops.GetMolFrags(mol, asMols=True))
        )
        if len(fragments) <= 1:
            return fragments[0] if fragments else None

        resolved_mode = mode or ("largest" if keep_largest else "reject")
        if resolved_mode == "keep":
            return mol
        if resolved_mode == "reject":
            logger.warning(
                f"{Chem.MolToSmiles(mol)} contains multiple unresolved fragments. Molecule rejected.")
            return None
        if resolved_mode != "largest":
            raise ValueError("Mixture mode must be 'keep', 'reject', or 'largest'.")

        candidates = [fragment for fragment in fragments if fragment.GetNumHeavyAtoms() > hac_threshold]
        if not candidates:
            logger.warning(
                f"{Chem.MolToSmiles(mol)} contains no fragments with >" + str(hac_threshold) + " heavy atoms")
            return None
        largest_size = max(fragment.GetNumHeavyAtoms() for fragment in candidates)
        largest = {
            ChemistryProcessor.standardize_smiles(fragment): fragment
            for fragment in candidates if fragment.GetNumHeavyAtoms() == largest_size
        }
        if len(largest) > 1:
            raise AmbiguousFragmentError(largest_size, list(largest))
        return next(iter(largest.values()))

    @staticmethod
    def remove_inorganic(mol: Chem.Mol) -> Optional[Chem.Mol]:
        """Reject carbon-free and selected small inorganic structures."""
        has_carbon = any(atom.GetSymbol() == 'C' for atom in mol.GetAtoms())
        if not has_carbon:
            return None

        inorganic_smarts = [
            '[#6]#[#7]', '[#8]=[#6]=[#8]', '[#6]#[#8]', '[#8]=[#6]=[#16]', '[#8]=[#6]=[#6]=[#6]=[#8]',
            '[#8]~[#6](=[#8])~[#8]', '[#16]=[#6]=[#16]', '[#7]=[#6]=[#8]', '[#8]-[#6]#[#7]',
            '[#16]-[#6]#[#7]', '[#7]=[#6]=[#16]', '[#34]-[#6]#[#7]', '[F,Cl,Br,I]-[#6]#[#7]',
            '[#7]#[#6]-[#6]#[#7]', '[Cl]-[#6](=[#8])-[Cl]', '[Cl]-[#6](=[#16])-[Cl]',
        ]
        for smarts in inorganic_smarts:
            pattern = Chem.MolFromSmarts(smarts)
            if not pattern:
                continue
            if mol.GetNumAtoms() == pattern.GetNumAtoms() and mol.HasSubstructMatch(pattern):
                return None

        return mol

    @staticmethod
    def reject_metals(mol: Chem.Mol) -> Optional[Chem.Mol]:
        """Reject structures containing a metal while leaving atom validation independent."""
        return None if ChemistryProcessor._contains_metal(mol) else mol

    def effective_atom(self, mol: Chem.Mol) -> Optional[Chem.Mol]:
        """Reject the whole molecule if any element is outside the allowed list."""
        invalid_elements = {atom.GetSymbol() for atom in mol.GetAtoms()} - self._valid
        if not invalid_elements:
            return mol
        logger.warning(
            f"Invalid elements {sorted(invalid_elements)} detected in molecule "
            f"{Chem.MolToSmiles(mol)}. Molecule rejected."
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
        """Get the same conservative single-organic pipeline used by DiptoxPipeline."""
        return cls.create_pipeline(
            processor_instance.reject_dummy_atoms,
            processor_instance.remove_isotopes,
            processor_instance.remove_hydrogens,
            lambda mol: processor_instance.remove_salts(mol),
            processor_instance.remove_solvents,
            lambda mol: processor_instance.remove_mixtures(mol, keep_largest=False),
            processor_instance.remove_inorganic,
            processor_instance.reject_radicals,
            processor_instance.neutralize_charges,
            processor_instance.collapse_identical_components,
        )
