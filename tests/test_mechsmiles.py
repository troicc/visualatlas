import csv
from pathlib import Path

import pytest

from visualatlas.mechsmiles import Arrow, canonical, parse, parse_arrow

FIXTURES = list(csv.DictReader(open(Path(__file__).parent / "data" / "uspto31k_steps.csv")))


def test_arrow_kinds():
    assert parse_arrow("(1, 2)") == Arrow("lp", (1,), 2)
    assert parse_arrow("((2,3), 3)") == Arrow("bond", (2, 3), 3)
    with pytest.raises(ValueError):
        parse_arrow("(1, 2, 3)")


def test_parse_ammonia_acetone_example():
    s = parse("C[C:2](=[O:3])C.[NH3:1]|(1,2);((2,3),3)")
    assert s.arrows == [Arrow("lp", (1,), 2), Arrow("bond", (2, 3), 3)]
    assert s.mol.GetAtomWithIdx(s.idx(1)).GetSymbol() == "N"
    assert s.mol.GetAtomWithIdx(s.idx(1)).GetTotalNumHs() == 3


def test_arrow_on_unmapped_atom_is_rejected():
    with pytest.raises(KeyError):
        parse("C[C:2](=O)C.[NH3:1]|(1,5)")


def test_state_without_arrows():
    assert parse("CCO").arrows == []


@pytest.mark.parametrize("row", FIXTURES, ids=[r["rxn_idx"] for r in FIXTURES])
def test_all_fixture_steps_parse(row):
    s = parse(row["step"])
    assert s.arrows
    assert canonical(row["expected_product"]) == row["expected_product"]


def test_fixture_products_differ_from_chrimp_only_by_stereo():
    from rdkit import Chem

    def strip(smi):
        m = Chem.MolFromSmiles(smi)
        Chem.RemoveStereochemistry(m)
        return Chem.MolToSmiles(m)

    dropped = 0
    for r in FIXTURES:
        assert strip(r["expected_product"]) == strip(r["chrimp_product"])
        dropped += int(r["stereo_dropped_by_chrimp"])
        assert (r["expected_product"] != r["chrimp_product"]) == bool(int(r["stereo_dropped_by_chrimp"]))
    assert dropped == 7


def test_canonical_ignores_atom_maps():
    assert canonical("[CH3:1][OH:2]") == canonical("CO") == "CO"
