"""Formal checker for elementary steps: arrows as electron moves on the molecular graph.

【标注】本模块由 AI（Claude）编写，2026-11-07。作者尚未逐行精读，面试前必须读懂并能复述每个设计取舍。
精读指南：docs/STUDY_checker.md。

Model: an electron ledger.
  - every atom a holds n(a) non-bonding electrons (lone pairs and unpaired electrons);
  - every bond a–b holds E(a,b) electrons (2 per bond order); bonds to implicit H stay fixed.
  - an arrow (visualatlas.moves.Move) takes 1 or 2 electrons from an atom or a bond and gives them to an
    atom or a bond. Nothing else changes; formal charge, bond order and unpaired electrons are derived:

        q(a) = V(a) - n(a) - sum_b o(a,b)          o(a,b) = E(a,b) / 2
        e(a) = n(a) + sum_b E(a,b) + 2 h(a)  =  V(a) - q(a) + sum_b o(a,b)   (electrons around a)
        unpaired(a) = n(a) mod 2

    V = valence electrons of the neutral atom (main group: group number mod 10; d block: group number),
    h = number of H atoms not drawn as separate atoms.

Because charge and bond order are derived from where the electrons are, conservation of atoms, electrons
and total charge holds by construction; what can go wrong is an arrow that takes electrons from where there
are none, an atom that ends with too many electrons, or a bond left with an odd number of electrons.

Rules:
  - 2-electron move from an atom: n(a) >= 2.          1-electron move from an atom: n(a) odd (it has an unpaired electron).
  - move from a bond: E(a,b) >= number of electrons moved.
  - after all arrows of the step (not after each arrow): every bond has an even number of electrons
    (two fish-hooks must meet to make a bond), and every atom satisfies e(a) <= limit:
        H, He 2; Li–Ne 8; main group period >= 3: 12; d-block metals: 18.
  - stereocentres whose bonds are not touched keep their configuration; stereo on touched atoms is not
    predicted (arrows do not say which face is attacked) and is not compared.
  - a drawn next state must have the same atoms and total charge as the state before, and must equal the
    computed one (canonical SMILES without atom maps); if it differs only in stereo, that is reported separately.
"""

from collections import Counter
from typing import NamedTuple

from rdkit import Chem

from ..mechsmiles import canonical, read_mol
from ..moves import moves_of

# ---------- issue codes ----------
NO_ELECTRONS = "no_electrons"  # an arrow takes electrons from an atom that has none (of the right kind)
NOT_BONDED = "not_bonded"  # an arrow takes electrons from a bond that does not hold enough
OVER_LIMIT = "over_limit"  # an atom ends the step with more electrons than its limit
HALF_BOND = "half_bond"  # a bond ends the step with an odd number of electrons
UNSUPPORTED = "unsupported"  # element without an electron limit rule (lanthanides, actinides, …)
ATOM_BALANCE = "atom_balance"  # drawn next state has a different atom multiset
CHARGE_BALANCE = "charge_balance"  # drawn next state has a different total charge
NEXT_MISMATCH = "next_mismatch"  # drawn next state is not what the arrows produce
STEREO_MISMATCH = "stereo_mismatch"  # same constitution, but an untouched stereocentre differs


class Issue(NamedTuple):
    code: str
    atoms: tuple  # atom-map numbers involved
    message: str


PT = Chem.GetPeriodicTable()
_D_BLOCK = {
    **{z: z - 18 for z in range(21, 31)},
    **{z: z - 36 for z in range(39, 49)},
    **{z: z - 68 for z in range(72, 81)},
}
_ORDER = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE, 4: Chem.BondType.QUADRUPLE}


# ---------- per-element constants ----------
def valence_electrons(atom):
    z = atom.GetAtomicNum()
    return _D_BLOCK.get(z, PT.GetNOuterElecs(z))


def electron_limit(atom):
    """Maximum electrons around the atom, or None when there is no rule for this element."""
    z = atom.GetAtomicNum()
    if z <= 2:
        return 2
    if z <= 10:
        return 8
    if z in _D_BLOCK:
        return 18
    if 57 <= z <= 71 or z >= 89:
        return None
    return 12


def _order_sum(mol, idx):
    """Sum of bond orders at atom idx, counting implicit/explicit H that are not separate atoms."""
    a = mol.GetAtomWithIdx(idx)
    return int(round(sum(b.GetBondTypeAsDouble() for b in a.GetBonds()))) + a.GetTotalNumHs()


def electron_count(mol, idx):
    """e(a) = V - q + sum of bond orders: electrons around the atom, shared ones counted fully."""
    a = mol.GetAtomWithIdx(idx)
    return valence_electrons(a) - a.GetFormalCharge() + _order_sum(mol, idx)


# ---------- the ledger ----------
class _Ledger:
    def __init__(self, step):
        self.step = step
        self.mol = Chem.RWMol(step.mol)
        for a in self.mol.GetAtoms():  # freeze H counts: RDKit must not re-derive them after bonds change
            a.SetNumExplicitHs(a.GetTotalNumHs())
            a.SetNoImplicit(True)
        self.n = {}  # atom idx -> non-bonding electrons
        for a in self.mol.GetAtoms():
            i = a.GetIdx()
            self.n[i] = valence_electrons(a) - a.GetFormalCharge() - _order_sum(self.mol, i)
        self.E = {}  # frozenset({i, j}) -> bond electrons
        for b in self.mol.GetBonds():
            self.E[frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))] = int(round(2 * b.GetBondTypeAsDouble()))
        self.E0 = dict(self.E)
        self.touched = set()
        self.issues = []

    def _ix(self, site):
        return tuple(self.step.idx(m) for m in site)

    def apply(self, mv):
        k = mv.electrons
        src, dst = self._ix(mv.src), self._ix(mv.dst)
        self.touched.update(src + dst)
        if len(src) == 1:
            (a,) = src
            # A fish-hook may also take one electron out of a lone pair (single-electron transfer, e.g. the
            # [O-][N+]=O <-> [O]N=O resonance in RMechDB), so the only requirement is enough electrons.
            if self.n[a] < k:
                what = "a lone pair" if k == 2 else "an electron"
                self.issues.append(Issue(NO_ELECTRONS, mv.src, f"atom {mv.src[0]} has no {what} to give"))
                return
            self.n[a] -= k
        else:
            key = frozenset(src)
            if self.E.get(key, 0) < k:
                self.issues.append(Issue(NOT_BONDED, mv.src, f"bond {mv.src} holds fewer than {k} electron(s)"))
                return
            self.E[key] -= k
        if len(dst) == 1:
            self.n[dst[0]] += k
        else:
            key = frozenset(dst)
            self.E[key] = self.E.get(key, 0) + k

    def finish(self):
        """Write the ledger back into a molecule; returns the product Mol (atom maps kept)."""
        mol = self.mol
        for key, e in self.E.items():
            i, j = tuple(key)
            order = e // 2
            b = mol.GetBondBetweenAtoms(i, j)
            if self.E0.get(key) == e:
                continue
            if b is not None and order == 0:
                mol.RemoveBond(i, j)
            elif b is not None:
                b.SetBondType(_ORDER[min(order, 4)])  # SetBondType keeps neighbour order, so chiral tags stay valid
                b.SetStereo(Chem.BondStereo.STEREONONE)
            elif order > 0:
                mol.AddBond(i, j, _ORDER[min(order, 4)])
        for i, n in self.n.items():
            a = mol.GetAtomWithIdx(i)
            osum = sum(e // 2 for key, e in self.E.items() if i in key) + a.GetNumExplicitHs()
            a.SetFormalCharge(valence_electrons(a) - n - osum)
            if i in self.touched:
                a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        _clear_bond_stereo(mol, self.touched)
        out = mol.GetMol()
        out.UpdatePropertyCache(strict=False)
        # Unpaired electrons follow RDKit's convention, the same one used when a drawn state is parsed:
        # electrons beyond the default valence count as radicals ([CH3] 1, [C]=O 2, [O] 1).
        Chem.AssignRadicals(out)
        return out


def _clear_bond_stereo(mol, atoms):
    """E/Z of a double bond with an end atom the arrows touched is not defined by the arrows.

    Only such double bonds lose their stereo, together with the '/' '\\' marks on their own neighbouring
    single bonds; an untouched double bond keeps its marks even if a touched atom is one of its neighbours.
    """
    for b in mol.GetBonds():
        if b.GetBondType() != Chem.BondType.DOUBLE:
            continue
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i not in atoms and j not in atoms:
            continue
        b.SetStereo(Chem.BondStereo.STEREONONE)
        for end in (i, j):
            for nb in mol.GetAtomWithIdx(end).GetBonds():
                if nb.GetIdx() != b.GetIdx():
                    nb.SetBondDir(Chem.BondDir.NONE)


def _run(step):
    led = _Ledger(step)
    for mv in moves_of(step):
        led.apply(mv)
    return led, led.finish()


def apply_arrows(step):
    """Apply the step's arrows in order and return the new molecule (atom maps kept).

    Impossible arrows are skipped here; check_step reports them.
    """
    return _run(step)[1]


# ---------- comparison helpers ----------
def _formula(mol):
    c = Counter()
    for a in mol.GetAtoms():
        c[a.GetSymbol()] += 1
        c["H"] += a.GetTotalNumHs()
    return c


def _charge(mol):
    return sum(a.GetFormalCharge() for a in mol.GetAtoms())


def _key(mol, stereo=True):
    m = Chem.Mol(mol)
    if not stereo:
        Chem.RemoveStereochemistry(m)
    k = canonical(m)
    if k is None:  # RDKit cannot sanitise (e.g. unusual valence): fall back to the kekulised, unsanitised SMILES
        for a in m.GetAtoms():
            a.SetAtomMapNum(0)
        k = Chem.MolToSmiles(m)
    return k


def check_step(step, drawn_next=None):
    """All issues of one elementary step; an empty list means the step is formally valid.

    step: visualatlas.mechsmiles.Step (arrows may be MechSMILES Arrows or Moves).
    drawn_next: SMILES string or RDKit Mol of the next drawn state, or None.
    """
    led, prod = _run(step)
    issues = list(led.issues)
    rev = {i: step.mol.GetAtomWithIdx(i).GetAtomMapNum() for i in range(step.mol.GetNumAtoms())}
    for key, e in led.E.items():
        if e % 2:
            i, j = sorted(key)
            issues.append(Issue(HALF_BOND, (rev[i], rev[j]), f"bond {rev[i]}–{rev[j]} ends with {e} electrons"))
    for a in prod.GetAtoms():
        lim = electron_limit(a)
        if lim is None:
            issues.append(Issue(UNSUPPORTED, (a.GetAtomMapNum(),), f"no electron limit rule for {a.GetSymbol()}"))
            continue
        e = electron_count(prod, a.GetIdx())
        if e > lim:
            issues.append(Issue(OVER_LIMIT, (a.GetAtomMapNum(),), f"{a.GetSymbol()} has {e} electrons (limit {lim})"))
    if drawn_next is None:
        return issues
    nxt = read_mol(drawn_next) if isinstance(drawn_next, str) else drawn_next
    if _formula(nxt) != _formula(step.mol):
        issues.append(
            Issue(
                ATOM_BALANCE,
                (),
                f"atoms differ: {dict(_formula(step.mol) - _formula(nxt))} lost, "
                f"{dict(_formula(nxt) - _formula(step.mol))} gained",
            )
        )
    if _charge(nxt) != _charge(step.mol):
        issues.append(Issue(CHARGE_BALANCE, (), f"total charge {_charge(step.mol)} -> {_charge(nxt)}"))
    if _key(prod, stereo=False) != _key(nxt, stereo=False):
        issues.append(Issue(NEXT_MISMATCH, (), f"arrows give {_key(prod)}, drawn {_key(nxt)}"))
    else:
        # stereo at atoms the arrows touched is not predicted, so it is ignored on both sides
        touched_maps = {step.mol.GetAtomWithIdx(i).GetAtomMapNum() for i in led.touched} - {0}
        nxt_cmp = Chem.Mol(nxt)
        mapped = [a for a in nxt_cmp.GetAtoms() if a.GetAtomMapNum() in touched_maps]
        if mapped:  # the drawn state carries atom maps: blank its stereo at the touched atoms
            idx = {a.GetIdx() for a in mapped}
            for i in idx:
                nxt_cmp.GetAtomWithIdx(i).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
            _clear_bond_stereo(nxt_cmp, idx)
            comparable = True
        else:  # no maps: only compare when no stereocentre was touched
            comparable = not any(
                step.mol.GetAtomWithIdx(i).GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for i in led.touched
            )
        if comparable and _key(prod) != _key(nxt_cmp):
            issues.append(
                Issue(STEREO_MISMATCH, (), f"untouched stereo differs: arrows give {_key(prod)}, drawn {_key(nxt)}")
            )
    return issues
