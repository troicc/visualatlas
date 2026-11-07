"""Formal checker for elementary steps: arrows as edits to the molecular graph, with electron bookkeeping.

STATUS: interface only. The implementation is to be written by the author (see PLAN.md, week 2).
tests/test_check.py holds the expected behaviour; those tests are skipped until the functions below
stop raising NotImplementedError.

Electron count of atom a (shared electrons counted fully, as in the octet rule):

    e(a) = V(a) - q(a) + sum_b o(a, b)

V = valence electrons of the neutral atom (group number), q = formal charge, o = bond order.
Non-bonding electrons are n(a) = V(a) - q(a) - sum_b o(a, b); implicit H count toward the bond sum.

Proposed rules (change them if you disagree, and record why in DECISIONS.md):
  - 'lp' arrow (a,) -> b: needs n(a) >= 2; then o(a,b) += 1, q(a) += 1, q(b) -= 1.
  - 'bond' arrow (a, b) -> c: needs o(a,b) >= 1; then o(a,b) -= 1 and
        c == b: q(a) += 1, q(b) -= 1        c == a: q(b) += 1, q(a) -= 1
        else:   o(b,c) += 1, q(a) += 1, q(c) -= 1
    (arrows are applied in the order given; a bond whose order reaches 0 is removed)
  - After all arrows, every atom must satisfy e(a) <= limit:
        H: 2; period 2 (Li–Ne): 8; main group period >= 3: 12; transition metals: 18.
  - Atom multiset and total charge are conserved between the state before and the drawn next state.
  - If a drawn next state is given, it must equal the state obtained by applying the arrows
    (compared as canonical SMILES without atom maps).
"""

from typing import NamedTuple

# issue codes
NO_ELECTRONS = "no_electrons"  # a lone-pair arrow starts on an atom with no lone pair
NOT_BONDED = "not_bonded"  # a bond arrow starts on a bond that does not exist
OVER_LIMIT = "over_limit"  # an atom exceeds its electron limit after the step
ATOM_BALANCE = "atom_balance"  # drawn next state has a different atom multiset
CHARGE_BALANCE = "charge_balance"  # drawn next state has a different total charge
NEXT_MISMATCH = "next_mismatch"  # drawn next state is not what the arrows produce


class Issue(NamedTuple):
    code: str
    atoms: tuple  # atom-map numbers involved
    message: str


def electron_count(mol, idx):
    """e(a) for atom index idx (see module docstring)."""
    raise NotImplementedError


def apply_arrows(step):
    """Apply step.arrows to step.mol in order and return the new molecule (atom maps kept).

    Raises nothing for chemically impossible arrows: report those through check_step instead.
    """
    raise NotImplementedError


def check_step(step, drawn_next=None):
    """All issues of one elementary step; an empty list means the step is formally valid.

    step: visualatlas.mechsmiles.Step; drawn_next: SMILES string or RDKit Mol of the next drawn state, or None.
    """
    raise NotImplementedError
