"""The visualatlas mechanism file: one JSON format for the checker, the renderer and LLMs.

【标注】本文件由 AI（Claude）编写，作者待精读。

A mechanism is a list of steps. Each step gives the state *before* it (a SMILES with atom maps on every atom
an arrow touches) and its arrows as electron moves. The state after step k is the species of step k+1, or
"final" after the last step; when "final" is omitted it is computed from the arrows.

    {
      "format": "visualatlas-mechanism/0.1",
      "title": "SN2 of ethyl bromide",
      "steps": [
        {"species": "C[CH2:1][Br:2].[OH-:3]",
         "moves": [{"e": 2, "from": [3], "to": [3, 1]},
                   {"e": 2, "from": [1, 2], "to": [2]}],
         "label": "SN2"}
      ],
      "final": "CCO.[Br-]"
    }

A move takes e electrons (2: curved arrow, 1: fish-hook) from an atom [a] or a bond [a, b] and gives them to an
atom [c] or a bond [c, d]. Optional drawing hints per move: "side" (up/down/left/right, which way the arrow
bulges). Optional per step: "label", "conditions".
"""

import json
from dataclasses import dataclass, field

import jsonschema

from .check import apply_arrows, check_step
from .mechsmiles import Step, canonical, parse, read_mol
from .moves import Move, moves_of

FORMAT = "visualatlas-mechanism/0.1"

_SITE = {"type": "array", "items": {"type": "integer", "minimum": 1}, "minItems": 1, "maxItems": 2}
SCHEMA = {
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "title": "visualatlas mechanism",
    "type": "object",
    "required": ["steps"],
    "additionalProperties": False,
    "properties": {
        "format": {"const": FORMAT},
        "title": {"type": "string"},
        "source": {"type": "object", "description": "where the mechanism comes from (dataset, row, DOI)"},
        "steps": {
            "type": "array",
            "minItems": 1,
            "items": {
                "type": "object",
                "required": ["species", "moves"],
                "additionalProperties": False,
                "properties": {
                    "species": {"type": "string", "description": "state before the step; SMILES with atom maps"},
                    "moves": {
                        "type": "array",
                        "minItems": 1,
                        "items": {
                            "type": "object",
                            "required": ["e", "from", "to"],
                            "additionalProperties": False,
                            "properties": {
                                "e": {"enum": [1, 2], "description": "2 = curved arrow, 1 = fish-hook"},
                                "from": {**_SITE, "description": "[atom] or [atom, atom] (a bond)"},
                                "to": {**_SITE, "description": "[atom] or [atom, atom] (a bond)"},
                                "side": {"enum": ["up", "down", "left", "right"]},
                            },
                        },
                    },
                    "label": {"type": "string"},
                    "conditions": {"type": "string"},
                },
            },
        },
        "final": {"type": "string", "description": "state after the last step"},
    },
}


@dataclass
class MechStep:
    species: str
    moves: list
    hints: list = field(default_factory=list)  # one dict per move (may be empty)
    label: str = ""
    conditions: str = ""

    def step(self):
        """The checker's Step for this step."""
        return Step(self.species, read_mol(self.species), list(self.moves))


@dataclass
class Mechanism:
    steps: list
    title: str = ""
    final: str | None = None
    source: dict = field(default_factory=dict)

    # ---------- states ----------
    def states(self):
        """SMILES of every drawn state: before each step, then after the last one (computed if not given)."""
        out = [s.species for s in self.steps]
        out.append(self.final if self.final is not None else _smiles(apply_arrows(self.steps[-1].step())))
        return out

    # ---------- JSON ----------
    def to_dict(self):
        d = {"format": FORMAT}
        if self.title:
            d["title"] = self.title
        if self.source:
            d["source"] = self.source
        d["steps"] = []
        for s in self.steps:
            moves = []
            for i, mv in enumerate(s.moves):
                m = {"e": mv.electrons, "from": list(mv.src), "to": list(mv.dst)}
                if i < len(s.hints) and s.hints[i]:
                    m.update(s.hints[i])
                moves.append(m)
            sd = {"species": s.species, "moves": moves}
            if s.label:
                sd["label"] = s.label
            if s.conditions:
                sd["conditions"] = s.conditions
            d["steps"].append(sd)
        if self.final is not None:
            d["final"] = self.final
        return d

    def dumps(self, **kw):
        return json.dumps(self.to_dict(), ensure_ascii=False, **kw)


def _smiles(mol):
    from rdkit import Chem

    return Chem.MolToSmiles(mol)


def from_dict(d):
    """Validate against SCHEMA and build a Mechanism. Raises jsonschema.ValidationError on bad input."""
    jsonschema.validate(d, SCHEMA)
    steps = []
    for sd in d["steps"]:
        moves, hints = [], []
        for m in sd["moves"]:
            moves.append(Move(m["e"], tuple(m["from"]), tuple(m["to"])))
            hints.append({k: v for k, v in m.items() if k not in ("e", "from", "to")})
        steps.append(MechStep(sd["species"], moves, hints, sd.get("label", ""), sd.get("conditions", "")))
    return Mechanism(steps, d.get("title", ""), d.get("final"), d.get("source", {}))


def loads(text):
    return from_dict(json.loads(text))


def load(path):
    with open(path, encoding="utf-8") as fh:
        return from_dict(json.load(fh))


# ---------- converters ----------
def from_mechsmiles(steps, title="", final=None):
    """A list of MechSMILES strings (one per step, in forward order) -> Mechanism."""
    out = []
    for text in steps:
        s = parse(text)
        out.append(MechStep(s.smiles, moves_of(s)))
    return Mechanism(out, title, final)


def from_rmechdb(record, title=""):
    """One RMechDB Record -> a one-step Mechanism (its written product becomes "final")."""
    s = record.step
    return Mechanism([MechStep(s.smiles, list(s.arrows), label=record.kind)], title, record.product)


# ---------- checking ----------
def check_mechanism(mech):
    """[(step index, issues)] for every step; each step is compared with the next drawn state."""
    states = mech.states()
    return [(k, check_step(s.step(), states[k + 1])) for k, s in enumerate(mech.steps)]


def same_state(a, b):
    return canonical(a) == canonical(b)
