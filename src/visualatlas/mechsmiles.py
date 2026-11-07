"""Parse MechSMILES (Neukomm, Jončev, Schwaller, arXiv 2512.05722) into an RDKit molecule and a list of arrows.

A MechSMILES step is "<mapped SMILES>|<arrow>;<arrow>;...". Arrows refer to atom-map numbers:

    (a, b)        the lone pair on atom a attacks atom b          -> Arrow('lp', (a,), b)
    ((a, b), c)   the electrons of bond a–b move towards atom c   -> Arrow('bond', (a, b), c)
                  c == b: the pair collapses onto b (heterolysis); c == a: onto a;
                  otherwise bond a–b breaks and bond b–c forms.

Example (ammonia adding to acetone): "C[C:2](=[O:3])C.[NH3:1]|(1,2);((2,3),3)".

The molecule is parsed without sanitisation and kekulised, because arrows act on localised bonds and
intermediates in mechanism datasets are not always valid RDKit structures (e.g. pentacoordinate Si).
"""

import ast
from dataclasses import dataclass, field

from rdkit import Chem


@dataclass(frozen=True)
class Arrow:
    kind: str  # 'lp' or 'bond'
    src: tuple  # (a,) for a lone pair, (a, b) for a bond; atom-map numbers
    dst: int  # atom-map number


@dataclass
class Step:
    smiles: str
    mol: Chem.Mol
    arrows: list = field(default_factory=list)

    def idx(self, map_number):
        """RDKit atom index of the atom with this map number."""
        for a in self.mol.GetAtoms():
            if a.GetAtomMapNum() == map_number:
                return a.GetIdx()
        raise KeyError(f"no atom with map number {map_number} in {self.smiles}")


def parse_arrow(text):
    t = ast.literal_eval(text.strip())
    if isinstance(t, tuple) and len(t) == 2 and all(isinstance(x, int) for x in t):
        return Arrow("lp", (t[0],), t[1])
    if (
        isinstance(t, tuple)
        and len(t) == 2
        and isinstance(t[0], tuple)
        and len(t[0]) == 2
        and all(isinstance(x, int) for x in t[0])
        and isinstance(t[1], int)
    ):
        return Arrow("bond", tuple(t[0]), t[1])
    raise ValueError(f"not a MechSMILES arrow: {text!r}")


def read_mol(smiles):
    """Parse a mapped SMILES without sanitisation, with implicit H counts computed and bonds kekulised."""
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        raise ValueError(f"unparsable SMILES: {smiles}")
    mol.UpdatePropertyCache(strict=False)
    Chem.Kekulize(mol, clearAromaticFlags=True)
    return mol


def parse(mech_smi):
    """MechSMILES string -> Step. A string without '|' is a plain state with no arrows."""
    smiles, _, arrows = mech_smi.partition("|")
    step = Step(smiles, read_mol(smiles), [parse_arrow(a) for a in arrows.split(";") if a.strip()])
    for ar in step.arrows:
        for m in (*ar.src, ar.dst):
            step.idx(m)  # raises if an arrow names an unmapped atom
    return step


def canonical(mol_or_smiles):
    """Canonical SMILES without atom maps, for comparing two states. Returns None if RDKit cannot sanitise it."""
    m = Chem.MolFromSmiles(mol_or_smiles, sanitize=False) if isinstance(mol_or_smiles, str) else Chem.Mol(mol_or_smiles)
    for a in m.GetAtoms():
        a.SetAtomMapNum(0)
    m = Chem.MolFromSmiles(Chem.MolToSmiles(m))
    return Chem.MolToSmiles(m) if m is not None else None
