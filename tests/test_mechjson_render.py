import csv
import json
import xml.etree.ElementTree as ET
from pathlib import Path

import jsonschema
import pytest

from visualatlas import cli, mechjson
from visualatlas.check import NEXT_MISMATCH
from visualatlas.moves import Move
from visualatlas.render import render
from visualatlas.rmechdb import parse_row

SN2 = "C[CH2:1][Br:2].[OH-:3]|(3, 1);((1, 2), 2)"
FIXTURES = list(csv.DictReader(open(Path(__file__).parent / "data" / "uspto31k_steps.csv")))
SCHEMA_FILE = Path(__file__).parents[1] / "schema" / "mechanism.schema.json"


def test_schema_file_is_current():
    assert json.loads(SCHEMA_FILE.read_text()) == mechjson.SCHEMA


def test_roundtrip_and_moves():
    m = mechjson.from_mechsmiles([SN2], "SN2", "CCO.[Br-]")
    assert m.steps[0].moves == [Move(2, (3,), (3, 1)), Move(2, (1, 2), (2,))]
    m2 = mechjson.loads(m.dumps())
    assert m2.to_dict() == m.to_dict()


@pytest.mark.parametrize(
    "bad",
    [
        {"steps": []},
        {"steps": [{"species": "C", "moves": [{"e": 3, "from": [1], "to": [2]}]}]},
        {"steps": [{"species": "C", "moves": [{"e": 2, "from": [1, 2, 3], "to": [2]}]}]},
        {"steps": [{"species": "C", "moves": [{"e": 2, "from": [1], "to": [2], "colour": "red"}]}]},
    ],
)
def test_schema_rejects(bad):
    with pytest.raises(jsonschema.ValidationError):
        mechjson.from_dict(bad)


def test_hints_survive_roundtrip():
    d = {"steps": [{"species": "C[CH2:1][Br:2].[OH-:3]", "moves": [{"e": 2, "from": [3], "to": [3, 1], "side": "up"}]}]}
    assert mechjson.from_dict(d).to_dict()["steps"][0]["moves"][0]["side"] == "up"


def test_check_mechanism():
    assert mechjson.check_mechanism(mechjson.from_mechsmiles([SN2], final="CCO.[Br-]")) == [(0, [])]
    issues = mechjson.check_mechanism(mechjson.from_mechsmiles([SN2], final="CC[O-].Br"))[0][1]
    assert {i.code for i in issues} == {NEXT_MISMATCH}


def test_final_state_computed_when_missing():
    m = mechjson.from_mechsmiles([SN2])
    assert mechjson.same_state(m.states()[-1], "CCO.[Br-]")


# ---------- rendering ----------
@pytest.mark.parametrize(
    "mech",
    [
        mechjson.from_mechsmiles([SN2], "SN2", "CCO.[Br-]"),
        mechjson.from_rmechdb(parse_row(["C[O:1][O:2]C>>C[O:1].C[O:2] 1,2-1;1,2-2", "", "", "homolyze"])),
        mechjson.from_rmechdb(
            parse_row(["[CH3:1].[CH2:2]=[CH:3]C>>[CH3:1][CH2:2][CH:3]C 1-1,2;2,3-1,2;2,3-3", "", "", "addition"])
        ),
    ],
    ids=["sn2", "homolysis", "radical-addition"],
)
def test_simple_steps_render_without_violations(mech, tmp_path):
    S = render(mech, str(tmp_path / "m.svg"))
    assert S.issues() == []
    assert len(S.__dict__["curls"]) == sum(len(s.moves) for s in mech.steps)
    assert ET.parse(tmp_path / "m.svg").getroot().tag.endswith("svg")


@pytest.mark.slow  # about 8 minutes: the curved-arrow search is slow on large molecules (see PLAN.md)
def test_every_fixture_step_renders():
    for r in FIXTURES:
        S = render(mechjson.from_mechsmiles([r["step"]]))
        assert len(S.__dict__["curls"]) == len(mechjson.from_mechsmiles([r["step"]]).steps[0].moves)


# ---------- command line ----------
def test_cli_convert_check_draw(tmp_path, capsys):
    src = tmp_path / "steps.txt"
    src.write_text(SN2 + "\n")
    js = tmp_path / "m.json"
    assert cli.main(["convert", "mechsmiles", str(src), "-o", str(js)]) == 0
    assert cli.main(["check", str(js)]) == 0
    assert cli.main(["draw", str(js), str(tmp_path / "m.svg")]) == 0
    assert "0 geometric violation" in capsys.readouterr().out
