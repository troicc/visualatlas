"""Read RMechDB (Tavakoli et al., J. Chem. Inf. Model. 2023, 10.1021/acs.jcim.2c01359).

【标注】本文件由 AI（Claude）编写，作者待精读，见 docs/STUDY_checker.md。

Each row of the CSV files is: "<reactants>>><products> <arrow codes>", conditions, stage, step type.
Arrow codes are ';'-separated "source-sink"; source and sink are an atom map number "a" or a bond "a,b".
Every RMechDB arrow is a fish-hook: it moves one electron.

The data are CC BY-NC-ND: they are read from a local download and never copied into this repository.
"""

import csv
from dataclasses import dataclass

from .mechsmiles import Step, read_mol
from .moves import Move


def _site(text):
    return tuple(int(x) for x in text.split(","))


def parse_codes(codes):
    out = []
    for c in codes.split(";"):
        c = c.strip()
        if not c:
            continue
        src, dst = c.split("-")
        out.append(Move(1, _site(src), _site(dst)))
    return out


@dataclass
class Record:
    step: Step
    product: str  # product SMILES as written in RMechDB
    condition: str
    stage: str
    kind: str


def parse_row(row):
    rxn, codes = row[0].split(" ", 1)  # SMILES has no spaces; the codes may ("11-11,2; 2-11,2")
    reac, prod = rxn.split(">>")
    step = Step(reac, read_mol(reac), parse_codes(codes))
    for mv in step.arrows:
        for m in mv.src + mv.dst:
            step.idx(m)
    return Record(step, prod, *(row[1:4] + [""] * (3 - len(row[1:4]))))


def read(path):
    """Yield (line number, Record or the exception raised while parsing it)."""
    with open(path, newline="") as fh:
        for i, row in enumerate(csv.reader(fh), start=1):
            if not row:
                continue
            try:
                yield i, parse_row(row)
            except Exception as e:  # noqa: BLE001 - parsing failures are data to report, not crashes
                yield i, e
