"""Element classification for metal detection, independent of allowed-atom rules.

DiPTox uses an explicit periodic-family convention: alkali/alkaline-earth
metals, d-block metals (including group 12), lanthanides, actinides, and the
listed p-block metals. B, Si, Ge, As, Sb and Te are treated as metalloids;
halogens and noble gases, including Ts and Og, are outside the metal set.

This is a declared application policy, not a unique IUPAC classification.
Polonium is included: the RSC entry describes its classification ambiguity
and discusses classifying it as a metal from its electrical conductivity.
The RSC entries also classify elements 113--116 as metals; many of their
bulk properties remain unknown. Sources checked 2026-09-22:

* https://periodic-table.rsc.org/element/84/polonium
* https://periodic-table.rsc.org/element/32/germanium
* https://periodic-table.rsc.org/element/51/antimony
* https://periodic-table.rsc.org/element/52/tellurium
* https://periodic-table.rsc.org/element/113/nihonium
* https://periodic-table.rsc.org/element/114/flerovium
* https://periodic-table.rsc.org/element/115/moscovium
* https://periodic-table.rsc.org/element/116/livermorium

The removable alkali subset preserves the existing salt-removal policy;
it does not imply that every occurrence is automatically a counterion.
"""

from __future__ import annotations

from numbers import Integral


REMOVABLE_ALKALI_METALS: frozenset[int] = frozenset({3, 11, 19, 37, 55, 87})

METAL_ATOMIC_NUMBERS: frozenset[int] = frozenset().union(
    REMOVABLE_ALKALI_METALS,
    {4, 12, 20, 38, 56, 88},                # Alkaline-earth metals.
    range(21, 31), range(39, 49),           # Sc--Zn, Y--Cd.
    range(72, 81), range(104, 113),         # Hf--Hg, Rf--Cn.
    range(57, 72), range(89, 104),          # La--Lu, Ac--Lr.
    {13, 31, 49, 50, 81, 82, 83, 84},     # Al, Ga, In, Sn, Tl, Pb, Bi, Po.
    {113, 114, 115, 116},                  # Nh, Fl, Mc, Lv.
)


def is_metal(atomic_number: int) -> bool:
    """Whether an integer atomic number belongs to DiPTox's declared metal set.

    Dummy atoms (0), unknown atomic numbers, and noninteger inputs are not
    metals. This classification neither permits nor rejects allowed atoms.
    """
    return isinstance(atomic_number, Integral) and atomic_number in METAL_ATOMIC_NUMBERS
