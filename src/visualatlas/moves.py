"""Electron moves: the one representation every arrow format is converted into before checking.

【标注】本文件由 AI（Claude）编写，作者待精读，见 docs/STUDY_checker.md。

A move transfers 1 electron (fish-hook arrow) or 2 electrons (ordinary curved arrow)
from a source to a sink. Source and sink are each either one atom (a,) or one bond (a, b);
all indices are atom-map numbers.

    lone pair a attacks b          MechSMILES (a, b)          Move(2, (a,), (a, b))
    bond a–b collapses onto b      MechSMILES ((a, b), b)     Move(2, (a, b), (b,))
    bond a–b moves to form b–c     MechSMILES ((a, b), c)     Move(2, (a, b), (b, c))
    homolysis of a–b               RMechDB "a,b-a;a,b-b"      Move(1, (a, b), (a,)), Move(1, (a, b), (b,))
    radical recombination a + b    RMechDB "a-a,b;b-a,b"      Move(1, (a,), (a, b)), Move(1, (b,), (a, b))
"""

from dataclasses import dataclass

from .mechsmiles import Arrow


@dataclass(frozen=True)
class Move:
    electrons: int  # 1 or 2
    src: tuple  # (a,) atom or (a, b) bond
    dst: tuple  # (a,) atom or (a, b) bond

    def atoms(self):
        return tuple(dict.fromkeys(self.src + self.dst))


def from_arrow(arrow):
    """MechSMILES Arrow -> Move."""
    if isinstance(arrow, Move):
        return arrow
    if arrow.kind == "lp":
        (a,) = arrow.src
        return Move(2, (a,), (a, arrow.dst))
    if arrow.kind == "bond":
        a, b = arrow.src
        c = arrow.dst
        if c == b:
            return Move(2, (a, b), (b,))
        if c == a:
            return Move(2, (a, b), (a,))
        return Move(2, (a, b), (b, c))
    raise ValueError(f"unknown arrow kind {arrow.kind!r}")


def moves_of(step):
    return [from_arrow(a) if isinstance(a, Arrow) else a for a in step.arrows]
