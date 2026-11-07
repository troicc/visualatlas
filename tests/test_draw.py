import xml.etree.ElementTree as ET

from visualatlas import demo
from visualatlas.draw import Canvas, M, check, mol_issues
from visualatlas.fonts import adv


def acetone():
    m = M()
    m.a("c", 0, 0)
    m.go("me1", "c", 210)
    m.go("me2", "c", 330)
    m.go("o", "c", 90, "O", t=2)
    return m


def test_font_metrics_available():
    assert 0.5 < adv("H") < 0.9
    assert adv("i") < adv("W")


def test_clean_molecule_has_no_issues():
    m = acetone()
    hard, soft = mol_issues(m, "acetone")
    assert hard == [] and soft == []
    assert check(m, "acetone") == []


def test_bad_angle_is_reported():
    m = M()
    m.a("a", 0, 0)
    m.go("b", "a", 0)
    m.go("c", "b", 90)  # 90° in an open chain: should be 120°
    hard, _ = mol_issues(m, "bent")
    assert any("b" in x and "应 120°" in x for x in hard), hard


def test_bad_bond_length_is_reported():
    m = M()
    m.a("a", 0, 0)
    m.a("b", 1.4, 0)
    m.b("a", "b")
    hard, _ = mol_issues(m, "long")
    assert any("长 1.40" in x for x in hard), hard


def test_crossing_bonds_are_reported():
    m = M()
    m.a("a", 0, 0).a("b", 1, 1).a("c", 0, 1).a("d", 1, 0)
    m.b("a", "b").b("c", "d")
    assert any("交叉" in x for x in check(m, "x"))


def test_demo_step_is_clean_and_valid_svg(tmp_path):
    S = demo.build()
    assert S.issues() == []
    out = tmp_path / "demo.svg"
    S.save(str(out))
    root = ET.parse(out).getroot()
    assert root.tag.endswith("svg")
    assert len(S.__dict__["curls"]) == 2


def test_impossible_arrow_is_recorded_not_silently_drawn():
    # a lone pair drawn inside the label's own charge superscript: every arc passes through the label
    S = Canvas(360, 220, 40)
    S.mol("k", acetone(), 230, 130)
    nu = M()
    nu.a("o", 0, 0, "HO−", anchor=1)
    S.mol("nu", nu, 110, 170)
    S.lone("nu", "o", 0)
    S.push(("lp", "nu", "o"), ("atom", "k", "c"))
    assert any("无可行弧" in x for x in S.issues())
    assert len(S.__dict__["curls"]) == 1  # still drawn, so the failure is visible
