"""Canvas: place drawn molecules, curved arrows and reaction arrows on one SVG, and collect every check.

Conventions (ACS / Wiley journal figures): black lines on white; one accent colour for the electrons
that move; reaction arrows with a solid head 0.30 B long; equilibria as a pair of half-headed arrows.

Every molecule placed with Canvas.mol() is checked (bond lengths, angles, label overlaps, crossings,
wedges); every curved arrow drawn with Canvas.push() is placed by search under the hard checks in
arrows.py. Canvas.issues() returns everything that failed: an empty list means the figure has no
detected geometric violation.
"""

import math

from . import arrows as _arrows
from . import geomcheck as _G
from . import mol as M2
from .mol import K, adv, layout

INK = "#1A1A1A"
ACC = _arrows.ACC
GREY = "#6B6B6B"
LAT = "'Liberation Sans',Arial,sans-serif"


def f(v):
    return f"{v:.1f}"


def mol_issues(m, name):
    """Geometry issues of one drawn molecule: (hard, soft). Soft issues are on atoms marked persp / ring_ok."""
    atoms, bonds, persp = {}, [], set()
    tsr = {k for k, a in m.at.items() if a.get("ts_ring")}
    for k, a in m.at.items():
        if a.get("geom") or a.get("pi"):
            continue
        atoms[k] = (a["x"], a["y"], a.get("label"))
        if a.get("ring_ok") or a.get("persp"):
            persp.add(k)
    for e in m.bd:
        if e["a"] in atoms and e["b"] in atoms:
            t = e["t"]
            if e.get("free") and t not in _G.SKIP:
                t = "d"  # free-length bond: check the angle, not the length
            bonds.append((e["a"], e["b"], t))
    r = _G.check_generic(name, atoms, bonds, persp, tsr)
    # collinearity is exempt only for atoms marked persp (perspective faces, bridged rings), not for ring_ok
    persp_only = {k for k, a in m.at.items() if a.get("persp")}
    # metals count as neighbours for collinearity (Rh–N–C drawn straight is also wrong); π dummy atoms do not
    atoms_m = dict(atoms)
    bonds_m = list(bonds)
    for k, a in m.at.items():
        if a.get("geom"):
            atoms_m[k] = (a["x"], a["y"], a.get("label") or "M")
    for e in m.bd:
        if (e["a"] in atoms_m and e["b"] in atoms_m) and not (e["a"] in atoms and e["b"] in atoms):
            bonds_m.append((e["a"], e["b"], "d" if e["t"] in _G.SKIP else e["t"]))
    r += [x for x in _G.collinear(name, atoms_m, bonds_m, persp_only, min_n=2) if x not in r]
    hard = [
        x
        for x in r
        if ("共线" in x and not any(f": {k}(" in x for k in persp_only))
        or not any(f": {k}(" in x for k in persp)
        and not any(e.get("len_ok") and f"{e['a']}–{e['b']}" in x for e in m.bd)
    ]
    soft = [x for x in r if x not in hard]
    return hard, soft


class Canvas:
    def __init__(s, w, h, B=40):
        s.w, s.h, s.B = w, h, B
        s.parts = []
        s.mols = {}
        s.geom_hard, s.geom_soft, s.arrow_issues = [], [], []

    # ---------- molecules ----------
    def mol(s, name, m, x, y, colors=None, bcolors=None, ink=INK):
        """Place molecule m with its origin at pixel (x, y) and check its geometry."""
        hard, soft = mol_issues(m, name)
        # overlaps, crossings and wedge ends from mol.check (angles and lengths are already covered above)
        hard += [x for x in M2.check(m, name, s.B) if "夹角" not in x and " 长 " not in x]
        s.geom_hard += hard
        s.geom_soft += soft
        s.parts.append(M2.svg_mol(m, s.B, x, y, ink, colors, bcolors))
        s.mols[name] = (m, x, y)
        return s

    def place(s, name, m, cx, cy, **kw):
        """Place molecule m centred on pixel (cx, cy); returns its pixel bounding box."""
        x0, y0, x1, y1 = M2.bbox(m, s.B)
        ox = cx - (x0 + x1) / 2
        oy = cy - (y0 + y1) / 2
        s.mol(name, m, ox, oy, **kw)
        return (ox + x0, oy + y0, ox + x1, oy + y1)

    def P(s, name, k):
        m, x, y = s.mols[name]
        a = m.at[k]
        return (x + a["x"] * s.B, y + a["y"] * s.B)

    def box(s, name):
        m, x, y = s.mols[name]
        b = M2.bbox(m, s.B)
        return (x + b[0], y + b[1], x + b[2], y + b[3])

    # ---------- electrons ----------
    def lone(s, name, k, ang, col=ACC, sep=0.17, r=0.052):
        """Draw the lone pair that reacts as two dots just outside the atom's label, in direction ang
        (mathematical angle). A later push(('lp', name, k), ...) starts from these dots."""
        m, ox, oy = s.mols[name]
        a = m.at[k]
        B = s.B
        fs = K["fs"] * B
        x, y = ox + a["x"] * B, oy + a["y"] * B
        t = math.radians(ang)
        ux, uy = math.cos(t), -math.sin(t)
        if a.get("label"):
            g = [q for q in layout(a["label"], fs, a.get("anchor", 0)) if q[5] is not None and q[2] == 0]
            ch = g[a.get("anchor", 0)] if g and a.get("anchor", 0) < len(g) else (g[0] if g else None)
            w = (adv(ch[0]) * fs / 2) if ch else fs * 0.35  # half width of the anchor glyph only
            h = fs * 0.36
            d = abs(ux) * w + abs(uy) * h + 0.09 * B
        else:
            d = 0.22 * B
        cx, cy = x + ux * d, y + uy * d
        px, py = -uy, ux
        for sg in (1, -1):
            s.parts.append(
                f'<circle cx="{f(cx + px * sep * B / 2 * sg)}" cy="{f(cy + py * sep * B / 2 * sg)}" '
                f'r="{f(r * B)}" fill="{col}"/>'
            )
        s.__dict__.setdefault("lps", {})[(name, k)] = ((cx, cy), (ux, uy))
        return s

    def rad(s, name, k, ang, gap=0.36, col=INK):
        """Draw an unpaired electron next to atom k in direction ang; ('rad', name, k) arrows start here."""
        x, y = s.P(name, k)
        r = math.radians(ang)
        p = (x + math.cos(r) * gap * s.B, y - math.sin(r) * gap * s.B)
        s.parts.append(f'<circle cx="{f(p[0])}" cy="{f(p[1])}" r="{f(0.065 * s.B)}" fill="{col}"/>')
        s.__dict__.setdefault("rads", {})[(name, k)] = p
        return s

    def push(s, src, dst, **kw):
        """Curved arrow from src to dst (see arrows.py for the endpoint syntax)."""
        _arrows.push(s, src, dst, **kw)
        return s

    # ---------- text ----------
    def chem(s, x, y, txt, size=None, col=INK, anchor="middle", bold=False):
        """Chemical text with the label syntax of mol.layout (H2 → H₂, ^+ → superscript, ...)."""
        fs = size or K["fs"] * s.B
        g = layout(txt, fs, 0)
        x0 = min(q[1] for q in g)
        x1 = max(q[1] + adv(q[0]) * q[3] for q in g)
        off = {"middle": -(x0 + x1) / 2, "start": -x0, "end": -x1}[anchor]
        sp = "".join(
            f'<tspan x="{f(x + off + q[1])}" y="{f(y + q[2])}" font-size="{f(q[3])}"'
            f"{' font-style=' + chr(34) + 'italic' + chr(34) if q[4] else ''}>"
            f"{'&lt;' if q[0] == '<' else q[0]}</tspan>"
            for q in g
        )
        bw = ' font-weight="bold"' if bold else ""
        s.parts.append(f'<text fill="{col}" style="font-family:{LAT}"{bw}>{sp}</text>')
        s.__dict__.setdefault("tboxes", []).append((txt, (x + off + x0, y - fs * 0.75, x + off + x1, y + fs * 0.25)))
        return s

    # ---------- reaction arrows ----------
    def _head(s, x, y, ang, col, scale=1.0, half=0):
        B = s.B
        L = 0.30 * B * scale
        W = 0.10 * B * scale
        ux, uy = math.cos(ang), math.sin(ang)
        nx, ny = -uy, ux
        bx, by = x - ux * L, y - uy * L
        cx, cy = x - ux * L * 0.78, y - uy * L * 0.78
        if half == 0:
            return (
                f'<path d="M{f(x)} {f(y)} L{f(bx + nx * W)} {f(by + ny * W)} L{f(cx)} {f(cy)} '
                f'L{f(bx - nx * W)} {f(by - ny * W)}Z" fill="{col}"/>'
            )
        return f'<path d="M{f(x)} {f(y)} L{f(bx + nx * W * half)} {f(by + ny * W * half)} L{f(cx)} {f(cy)}Z" fill="{col}"/>'

    def arrow(s, x0, y0, x1, y1, col=INK):
        w = K["lw"] * s.B * 1.05
        a = math.atan2(y1 - y0, x1 - x0)
        L = 0.30 * s.B * 0.78
        s.parts.append(
            f'<line x1="{f(x0)}" y1="{f(y0)}" x2="{f(x1 - math.cos(a) * L)}" y2="{f(y1 - math.sin(a) * L)}" '
            f'stroke="{col}" stroke-width="{f(w)}"/>'
        )
        s.parts.append(s._head(x1, y1, a, col))
        s.__dict__.setdefault("apaths", []).append(
            (f"({x0:.0f},{y0:.0f})→({x1:.0f},{y1:.0f})", _samples((x0, y0), (x1, y1)))
        )
        return s

    def eq(s, x0, y, x1, col=INK):
        """Equilibrium arrow: upper line with a half head on the right, lower line with a half head on the left."""
        sep = 0.11 * s.B
        w = K["lw"] * s.B * 1.05
        L = 0.32 * s.B
        hw = 0.13 * s.B
        s.parts.append(
            f'<line x1="{f(x0)}" y1="{f(y - sep)}" x2="{f(x1)}" y2="{f(y - sep)}" stroke="{col}" stroke-width="{f(w)}"/>'
        )
        s.parts.append(
            f'<path d="M{f(x1 + w * 0.3)} {f(y - sep)} L{f(x1 - L)} {f(y - sep - hw)} '
            f'L{f(x1 - L * 0.8)} {f(y - sep)}Z" fill="{col}"/>'
        )
        s.parts.append(
            f'<line x1="{f(x0)}" y1="{f(y + sep)}" x2="{f(x1)}" y2="{f(y + sep)}" stroke="{col}" stroke-width="{f(w)}"/>'
        )
        s.parts.append(
            f'<path d="M{f(x0 - w * 0.3)} {f(y + sep)} L{f(x0 + L)} {f(y + sep + hw)} '
            f'L{f(x0 + L * 0.8)} {f(y + sep)}Z" fill="{col}"/>'
        )
        s.__dict__.setdefault("apaths", []).append(
            (f"⇌({x0:.0f},{y:.0f})", _samples((x0, y - sep), (x1, y - sep)) + _samples((x0, y + sep), (x1, y + sep)))
        )
        return s

    # ---------- output ----------
    def issues(s):
        """Every hard violation found so far (molecule geometry + curved arrows). Empty = clean figure."""
        return list(s.geom_hard) + list(s.arrow_issues)

    def svg(s):
        return (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{s.w}" height="{s.h}" viewBox="0 0 {s.w} {s.h}">'
            f'<rect width="{s.w}" height="{s.h}" fill="#FFFFFF"/>' + "\n".join(s.parts) + "</svg>"
        )

    def save(s, path):
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(s.svg())
        return path


def _samples(p0, p1, n=30, trim=6):
    L = math.dist(p0, p1) or 1
    ux, uy = (p1[0] - p0[0]) / L, (p1[1] - p0[1]) / L
    a = (p0[0] + ux * trim, p0[1] + uy * trim)
    b = (p1[0] - ux * trim, p1[1] - uy * trim)
    return [(a[0] + (b[0] - a[0]) * k / n, a[1] + (b[1] - a[1]) * k / n) for k in range(n + 1)]
