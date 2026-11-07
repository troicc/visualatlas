import pytest

from visualatlas.draw import M
from visualatlas.stereo import cip_from_drawing

ALA = "[NH2:1][CH:2]([CH3:3])[C:4](=O)O"
KEYMAP = {1: "n", 2: "ca", 3: "me", 4: "c"}


def alanine(methyl_bond):
    m = M()
    m.a("ca", 0, 0)
    m.go("n", "ca", 210, "H2N")
    m.go("c", "ca", 330)
    m.go("me", "ca", 90, t=methyl_bond)
    return m


@pytest.mark.parametrize("bond, expected", [("w", "S"), ("h", "R")])
def test_cip_read_back_from_wedges(bond, expected):
    # N lower left, CO2H lower right, CH3 up: N → CO2H → CH3 runs anticlockwise.
    # With CH3 wedged the implicit H points away from the viewer, so the centre is S (L-alanine).
    assert cip_from_drawing(alanine(bond), ALA, KEYMAP, "ca") == expected


def test_unmapped_neighbour_is_rejected():
    with pytest.raises(ValueError):
        cip_from_drawing(alanine("w"), ALA, {1: "n", 2: "ca", 3: "me"}, "ca")
