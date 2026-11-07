"""Fish-hook arrows, half bonds, stereo comparison and metal electron counts.

All examples are textbook steps written for these tests (RMechDB rows are CC BY-NC-ND and are not copied here).
Steps use RMechDB's notation: "<reactants>>><products> <source>-<sink>;...", one electron per arrow.
"""

import pytest

from visualatlas.check import (
    HALF_BOND,
    NO_ELECTRONS,
    STEREO_MISMATCH,
    apply_arrows,
    check_step,
    electron_count,
)
from visualatlas.mechsmiles import canonical, parse
from visualatlas.rmechdb import parse_codes, parse_row
from visualatlas.moves import Move


def rec(text):
    return parse_row([text, "", "", ""])


def test_codes_with_spaces_after_semicolons():
    assert parse_codes("11-11,2; 2-11,2") == [Move(1, (11,), (11, 2)), Move(1, (2,), (11, 2))]


@pytest.mark.parametrize(
    "text",
    [
        "C[O:1][O:2]C>>C[O:1].C[O:2] 1,2-1;1,2-2",  # peroxide homolysis
        "[CH3:1].[CH3:2]>>[CH3:1][CH3:2] 1-1,2;2-1,2",  # radical recombination
        "C[O:1].[H:2][CH2:3]C>>C[O:1][H:2].[CH2:3]C 1-1,2;2,3-1,2;2,3-3",  # H-atom abstraction
        "[CH3:1].[CH2:2]=[CH2:3]>>[CH3:1][CH2:2][CH2:3] 1-1,2;2,3-1,2;2,3-3",  # radical addition to an alkene
        "[CH2:1]=[CH:2][CH2:3]>>[CH2:1][CH:2]=[CH2:3] 3-2,3;1,2-2,3;1,2-1",  # allyl radical resonance
        "[O-:1][N+:2]=O>>[O:1][N:2]=O 1-2",  # NO2: one electron shifts between its two textbook resonance drawings
    ],
)
def test_valid_radical_steps(text):
    r = rec(text)
    assert check_step(r.step, r.product) == []
    assert canonical(apply_arrows(r.step)) == canonical(r.product)


def test_one_fish_hook_leaves_a_half_bond():
    r = rec("[CH3:1].[CH3:2]>>CC 1-1,2")
    assert HALF_BOND in {i.code for i in check_step(r.step, r.product)}


def test_fish_hook_from_atom_without_electrons():
    r = rec("C[O:1].[H:2][CH2:3]C>>CO.[CH2]C 1-1,2;2-1,2;2,3-3")  # should start from the H–C bond, not from H
    issues = check_step(r.step, r.product)
    assert any(i.code == NO_ELECTRONS and 2 in i.atoms for i in issues)


# ---------- stereo ----------
SN2_STEREO = "C[C@H](O)C[CH2:1][Br:2].[OH-:3]|(3, 1);((1, 2), 2)"


def test_untouched_stereocentre_kept():
    assert check_step(parse(SN2_STEREO), "C[C@H](O)CCO.[Br-]") == []


def test_untouched_stereocentre_lost_in_drawn_state():
    assert {i.code for i in check_step(parse(SN2_STEREO), "CC(O)CCO.[Br-]")} == {STEREO_MISMATCH}


def test_untouched_stereocentre_inverted_in_drawn_state():
    assert {i.code for i in check_step(parse(SN2_STEREO), "C[C@@H](O)CCO.[Br-]")} == {STEREO_MISMATCH}


def test_stereo_of_new_double_bond_is_not_judged():
    # the arrows cannot say E or Z; a drawn E isomer is accepted
    r = rec("[CH2:1]=[CH:2][CH:3]C>>[CH2:1]/[CH:2]=[CH:3]/C 3-2,3;1,2-2,3;1,2-1")
    assert check_step(r.step, r.product) == []


def test_untouched_double_bond_next_to_touched_atom_keeps_its_geometry():
    # 2-methylbut-2-enyl radical + O2; the drawn product marks the geometry on the other substituent
    r = rec("[CH2:2]/C(C)=C/CC.[O][O:1]>>C/C([CH2:2][O:1][O])=C\\CC 1-1,2;2-2,1")
    assert check_step(r.step, r.product) == []


# ---------- metals ----------
@pytest.mark.parametrize(
    "smiles, mapnum, expected",
    [
        ("[Cu:1]Cl", 1, 12),  # Cu(I): d10 + one bond
        # trans-PdCl2(PH3)2, a 16-electron complex. Dative P->Pd bonds must be written charge-separated
        # ([PH3+], [Pd-2]); written as plain single bonds the P atoms would be wrong, not the Pd count.
        ("Cl[Pd-2:1](Cl)([PH3+])[PH3+]", 1, 16),
        ("[Fe:1]", 1, 8),
    ],
)
def test_metal_electron_counts(smiles, mapnum, expected):
    s = parse(smiles)
    assert electron_count(s.mol, s.idx(mapnum)) == expected
