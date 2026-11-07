"""Demo: one hand-written elementary step (hydroxide adding to the carbonyl of acetone), checked and saved as SVG.

visualatlas-demo out.svg
"""

import sys

from .draw import Canvas, M


def acetone():
    m = M()
    m.a("c", 0, 0)
    m.go("me1", "c", 210)
    m.go("me2", "c", 330)
    m.go("o", "c", 90, "O", t=2)
    return m


def hydroxide():
    m = M()
    m.a("o", 0, 0, "HO−", anchor=1)
    return m


def build(B=40):
    S = Canvas(360, 220, B)
    S.mol("k", acetone(), 230, 130)
    S.mol("nu", hydroxide(), 105, 175)
    S.lone("nu", "o", -45)
    S.push(("lp", "nu", "o"), ("atom", "k", "c"))
    S.push(("bond", "k", "c", "o"), ("atom", "k", "o"))
    return S


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    out = argv[0] if argv else "demo.svg"
    S = build()
    S.save(out)
    problems = S.issues()
    for p in problems:
        print("  ✗", p)
    print(f"{out}: {len(problems)} issue(s)")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
