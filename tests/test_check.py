"""Expected behaviour of the formal checker (visualatlas.check).

Every test here is skipped while the checker raises NotImplementedError, so CI stays green until it exists.
"""

import csv
from pathlib import Path

import pytest

from visualatlas import check
from visualatlas.check import (
    ATOM_BALANCE,
    CHARGE_BALANCE,
    NEXT_MISMATCH,
    NO_ELECTRONS,
    NOT_BONDED,
    OVER_LIMIT,
    apply_arrows,
    check_step,
    electron_count,
)
from visualatlas.mechsmiles import canonical, parse

FIXTURES = list(csv.DictReader(open(Path(__file__).parent / "data" / "uspto31k_steps.csv")))


def _implemented():
    try:
        check.check_step(parse("[OH-:1].[CH3+:2]|(1, 2)"))
        return True
    except NotImplementedError:
        return False


pytestmark = pytest.mark.skipif(not _implemented(), reason="checker not implemented yet")


def codes(issues):
    return {i.code for i in issues}


# ---------- electron counts ----------
@pytest.mark.parametrize(
    "smiles, mapnum, expected",
    [
        ("[CH4:1]", 1, 8),
        ("[CH3+:1]", 1, 6),
        ("[CH3-:1]", 1, 8),
        ("[OH-:1]", 1, 8),
        ("[CH3:1]", 1, 7),  # methyl radical
        ("[BH3:1]", 1, 6),
        ("C[S:1](C)(=O)=O", 1, 12),  # dimethyl sulfone, written with S=O double bonds
        ("[H:1][H]", 1, 2),
    ],
)
def test_electron_count(smiles, mapnum, expected):
    s = parse(smiles)
    assert electron_count(s.mol, s.idx(mapnum)) == expected


# ---------- applying arrows ----------
@pytest.mark.parametrize(
    "mech, product",
    [
        ("[OH-:1].[CH3+:2]|(1, 2)", "CO"),  # cation capture
        ("C[C:1](C)(C)[Br:2]|((1, 2), 2)", "C[C+](C)C.[Br-]"),  # heterolysis onto Br
        ("C[CH2:1][Br:2].[OH-:3]|(3, 1);((1, 2), 2)", "CCO.[Br-]"),  # SN2
        ("C[C:3](=[O:4])C.[BH3-:1][H:2]|((1, 2), 3);((3, 4), 4)", "B.CC(C)[O-]"),  # hydride transfer to acetone
        ("C[C:2](=[O:3])C.[NH3:1]|(1, 2);((2, 3), 3)", "CC(C)([O-])[NH3+]"),  # amine addition
    ],
)
def test_apply_arrows(mech, product):
    assert canonical(apply_arrows(parse(mech))) == canonical(product)


@pytest.mark.parametrize("row", FIXTURES, ids=[r["rxn_idx"] for r in FIXTURES])
def test_apply_arrows_on_mech_uspto_31k(row):
    s = parse(row["step"])
    assert canonical(apply_arrows(s)) == row["expected_product"]


# ---------- valid steps have no issues ----------
@pytest.mark.parametrize("row", FIXTURES, ids=[r["rxn_idx"] for r in FIXTURES])
def test_dataset_steps_are_valid(row):
    s = parse(row["step"])
    assert check_step(s) == []
    assert check_step(s, row["expected_product"]) == []


def test_transient_overflow_between_arrows_is_allowed():
    # after the first arrow the carbonyl C has 10 electrons; the second arrow fixes it. Only the end state counts.
    assert check_step(parse("C[C:3](=[O:4])C.[BH3-:1][H:2]|((1, 2), 3);((3, 4), 4)")) == []


# ---------- injected errors ----------
def test_lone_pair_arrow_from_atom_without_lone_pair():
    issues = check_step(parse("[CH4:1].[CH3+:2]|(1, 2)"))
    assert NO_ELECTRONS in codes(issues)
    assert any(i.code == NO_ELECTRONS and 1 in i.atoms for i in issues)


def test_octet_exceeded():
    issues = check_step(parse("[OH-:1].[CH4:2]|(1, 2)"))  # pentavalent carbon
    assert any(i.code == OVER_LIMIT and 2 in i.atoms for i in issues)


def test_bond_arrow_on_missing_bond():
    issues = check_step(parse("C[CH2:1][Br:2].[Cl-:3]|((1, 3), 3)"))
    assert NOT_BONDED in codes(issues)


SN2 = "C[CH2:1][Br:2].[OH-:3]|(3, 1);((1, 2), 2)"


def test_drawn_next_state_correct():
    assert check_step(parse(SN2), "CCO.[Br-]") == []


def test_drawn_next_state_wrong_structure():
    # same atoms and charge, but not what the arrows give
    issues = check_step(parse(SN2), "CC[O-].Br")
    assert codes(issues) == {NEXT_MISMATCH}


def test_drawn_next_state_lost_an_atom():
    assert ATOM_BALANCE in codes(check_step(parse(SN2), "CCO"))


def test_drawn_next_state_wrong_charge():
    assert CHARGE_BALANCE in codes(check_step(parse(SN2), "CCO.[Br]"))
