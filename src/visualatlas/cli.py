"""Command line: convert to the mechanism JSON, check it, draw it.

【标注】本文件由 AI（Claude）编写，作者待精读。

    visualatlas convert mechsmiles steps.txt -o mech.json   # one MechSMILES per line, in forward order
    visualatlas check mech.json                              # exit code 1 if any step has an issue
    visualatlas draw mech.json out.svg                       # prints the geometric violations
"""

import argparse
import sys

from rdkit import RDLogger

from . import mechjson


def _convert(a):
    lines = [ln.strip() for ln in open(a.input, encoding="utf-8") if ln.strip() and not ln.startswith("#")]
    if a.kind == "mechsmiles":
        mech = mechjson.from_mechsmiles(lines, title=a.title or "")
    else:
        from . import rmechdb

        rec = rmechdb.parse_row([lines[0], "", "", ""])
        mech = mechjson.from_rmechdb(rec, title=a.title or "")
    text = mech.dumps(indent=2)
    if a.output:
        with open(a.output, "w", encoding="utf-8") as fh:
            fh.write(text + "\n")
    else:
        print(text)
    return 0


def _check(a):
    mech = mechjson.load(a.input)
    bad = 0
    for k, issues in mechjson.check_mechanism(mech):
        for i in issues:
            bad += 1
            print(f"step {k + 1}: {i.code} {i.atoms} {i.message}")
    print(f"{len(mech.steps)} step(s), {bad} issue(s)")
    return 1 if bad else 0


def _draw(a):
    from .render import render

    S = render(mechjson.load(a.input), a.output)
    for x in S.issues():
        print("  ✗", x)
    print(f"{a.output}: {len(S.issues())} geometric violation(s)")
    return 0


def main(argv=None):
    RDLogger.DisableLog("rdApp.*")
    p = argparse.ArgumentParser(prog="visualatlas")
    sub = p.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("convert", help="MechSMILES or RMechDB text -> mechanism JSON")
    c.add_argument("kind", choices=["mechsmiles", "rmechdb"])
    c.add_argument("input")
    c.add_argument("-o", "--output")
    c.add_argument("--title")
    c.set_defaults(fn=_convert)
    k = sub.add_parser("check", help="formal check of every step")
    k.add_argument("input")
    k.set_defaults(fn=_check)
    d = sub.add_parser("draw", help="render to SVG")
    d.add_argument("input")
    d.add_argument("output")
    d.set_defaults(fn=_draw)
    a = p.parse_args(argv)
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
