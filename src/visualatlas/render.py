"""Render a Mechanism (mechjson) as one SVG: every state in a row, the step's arrows drawn on the state before it.

【标注】本文件由 AI（Claude）编写，作者待精读。

Version 0 (the ablation baseline ② of PLAN.md, experiment 1): fragment coordinates come straight from RDKit
(CoordGen), fragments of one state are placed left to right, and only the curved arrows are placed by search
(draw/arrows.py). No layout search yet.

Drawing conventions (ACS 1996):
  - carbon is a bare vertex unless it is alone, charged or carries an unpaired electron next to nothing else;
  - heteroatoms are labelled with their H count, the H written on the side away from the bonds (OH / HO);
  - charges: superscript in the label, or a small sign next to a bare carbon vertex;
  - unpaired electrons: a dot on the freest side of the atom;
  - lone pairs are drawn only where an arrow starts from them.
"""

import math

from rdkit import Chem
from rdkit.Chem import rdDepictor

from .draw import Canvas, M
from .mechsmiles import read_mol

rdDepictor.SetPreferCoordGen(True)

B = 40  # pixels per bond length
GAP_FRAG = 1.6  # bond lengths between fragments of one state
GAP_PANEL = 3.6  # bond lengths for the reaction arrow between states
PAD = 1.4


# ---------- one state -> drawn fragments ----------
class Fragment:
    def __init__(s, name, m, keys):
        s.name, s.m, s.keys = name, m, keys  # keys: rdkit atom idx -> atom key in m


def _label(atom, right_heavy):
    """Label text for an atom, or None for a bare carbon vertex."""
    sym = atom.GetSymbol()
    h = atom.GetTotalNumHs()
    q = atom.GetFormalCharge()
    rad = atom.GetNumRadicalElectrons()
    if sym == "C" and atom.GetDegree() > 0 and (atom.GetDegree() > 1 or (q == 0 and rad == 0)):
        return None
    hs = "" if h == 0 else ("H" if h == 1 else f"H{h}")
    chg = "" if q == 0 else (("+" if q > 0 else "−") if abs(q) == 1 else f"{abs(q)}{'+' if q > 0 else '−'}")
    if chg and len(chg) > 1:
        chg = "^{" + chg + "}"
    core = (hs + sym) if (right_heavy and hs) else (sym + hs)
    return core + chg


def _free_angle(m, k, prefer=None):
    """Mathematical angle (deg) of the largest gap between the bonds at atom k (optionally biased to 'prefer')."""
    a = m.at[k]
    angs = []
    for e in m.bd:
        if k in (e["a"], e["b"]):
            o = m.at[e["b"] if e["a"] == k else e["a"]]
            angs.append(math.degrees(math.atan2(-(o["y"] - a["y"]), o["x"] - a["x"])) % 360)
    if not angs:
        return prefer if prefer is not None else 90.0
    best, score = 0.0, -1e9
    for th in range(0, 360, 10):
        gap = min(min(abs(th - g) % 360, 360 - abs(th - g) % 360) for g in angs)
        sc = gap + (25 * math.cos(math.radians(th - prefer)) if prefer is not None else 0)
        if sc > score:
            best, score = float(th), sc
    return best


def state_fragments(smiles, prefix):
    """Mapped SMILES -> list of Fragment with coordinates in bond lengths (y down), plus map number -> (name, key)."""
    mol = read_mol(smiles)
    Chem.GetSymmSSSR(mol)
    rdDepictor.Compute2DCoords(mol)
    conf = mol.GetConformer()
    # normalise to bond length 1 (RDKit's depictors differ: 1.5 classic, 1.0 CoordGen)
    lens = sorted(
        (conf.GetAtomPosition(b.GetBeginAtomIdx()) - conf.GetAtomPosition(b.GetEndAtomIdx())).Length()
        for b in mol.GetBonds()
    )
    unit = lens[len(lens) // 2] if lens else 1.0
    try:
        Chem.WedgeMolBonds(mol, conf)
    except Exception:  # noqa: BLE001 - unsanitised intermediates may not support wedging
        pass
    ri = mol.GetRingInfo()
    frags, where = [], {}
    for j, idxs in enumerate(Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False)):
        name = f"{prefix}f{j}"
        m = M()
        keys = {}
        for i in idxs:
            p = conf.GetAtomPosition(i)
            keys[i] = f"a{i}"
            atom = mol.GetAtomWithIdx(i)
            nb = [conf.GetAtomPosition(n.GetIdx()).x for n in atom.GetNeighbors()]
            right = bool(nb) and sum(nb) / len(nb) < p.x  # bonds come from the left -> write H after? no: H before
            lab = _label(atom, right_heavy=not right and bool(nb))
            m.a(keys[i], p.x / unit, -p.y / unit, lab)
            m.at[keys[i]]["_q"] = atom.GetFormalCharge()
            m.at[keys[i]]["_rad"] = atom.GetNumRadicalElectrons()
            if atom.GetAtomMapNum():
                where[atom.GetAtomMapNum()] = (name, keys[i])
        for b in mol.GetBonds():
            i, k = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            if i not in keys or k not in keys:
                continue
            order = int(round(b.GetBondTypeAsDouble())) or 1
            t = {1: 1, 2: 2, 3: 3}.get(order, 1)
            d = b.GetBondDir()
            if t == 1 and d == Chem.BondDir.BEGINWEDGE:
                t = "w"
            elif t == 1 and d == Chem.BondDir.BEGINDASH:
                t = "h"
            kw = {}
            if t == 2 and ri.NumBondRings(b.GetIdx()):
                ring = min((r for r in ri.AtomRings() if i in r and k in r), key=len)
                kw = dict(side="in", ring=[keys[x] for x in ring])
            m.b(keys[i], keys[k], t, **kw)
        frags.append(Fragment(name, m, keys))
    return frags, where


def _bbox(m):
    from .draw import bbox

    return bbox(m, 1)


# ---------- mechanism -> SVG ----------
def _arrow_for(move, where, frags_by_name):
    """Move -> (src, dst, kw) in canvas terms."""
    half = move.electrons == 1
    kw = {"half": half}

    def site(t):
        return tuple(x for m in t for x in where[m])  # (name, key) or (name, key, name, key)

    def bonded(t):
        (na, ka), (nb, kb) = where[t[0]], where[t[1]]
        if na != nb:
            return False
        m = frags_by_name[na].m
        return any({e["a"], e["b"]} == {ka, kb} for e in m.bd)

    # source
    if len(move.src) == 1:
        name, key = where[move.src[0]]
        atom = frags_by_name[name].m.at[key]
        src = ("rad", name, key) if (half and atom.get("_rad")) else ("lp", name, key)
    else:
        (na, ka), (_, kb) = where[move.src[0]], where[move.src[1]]
        src = ("bond", na, ka, kb)
    # sink
    if len(move.dst) == 1:
        dst = ("atom",) + site(move.dst)
    elif bonded(move.dst):
        (na, ka), (_, kb) = where[move.dst[0]], where[move.dst[1]]
        dst = ("bond", na, ka, kb)
    elif half:  # a new bond made of two fish-hooks: both meet between the atoms
        dst = ("mid",) + site(move.dst)
    else:  # a new bond made by a curved arrow: point at the atom that receives the bond
        other = [a for a in move.dst if a not in move.src] or [move.dst[-1]]
        dst = ("atom",) + where[other[0]]
    return src, dst, kw


def _stretch_label_bonds(m):
    """Bonds with a label at both ends are drawn 1.5 long (otherwise the labels eat the bond).
    A terminal atom is simply moved outwards; other cases are left for the layout search."""
    deg = {k: 0 for k in m.at}
    for e in m.bd:
        deg[e["a"]] += 1
        deg[e["b"]] += 1
    for e in m.bd:
        a, b = m.at[e["a"]], m.at[e["b"]]
        if not (a.get("label") and b.get("label")):
            continue
        for far, near in ((e["b"], e["a"]), (e["a"], e["b"])):
            if deg[far] == 1:
                fa, na = m.at[far], m.at[near]
                dx, dy = fa["x"] - na["x"], fa["y"] - na["y"]
                L = math.hypot(dx, dy) or 1
                fa["x"], fa["y"] = na["x"] + dx / L * 1.5, na["y"] + dy / L * 1.5
                break
    # label-label bonds inside a chain: relax the whole fragment with those bonds at 1.5, others at 1
    inner = [
        e
        for e in m.bd
        if m.at[e["a"]].get("label") and m.at[e["b"]].get("label") and deg[e["a"]] > 1 and deg[e["b"]] > 1
    ]
    if inner:
        from .draw import check, relax

        before = {k: (a["x"], a["y"]) for k, a in m.at.items()}
        n0 = len(check(m))
        for e in m.bd:
            e["L0"] = 1.5 if (m.at[e["a"]].get("label") and m.at[e["b"]].get("label")) else 1.0
        relax(m, anchor_w=0.01)
        for e in m.bd:
            e.pop("L0", None)
        if len(check(m)) > n0:  # relaxation made it worse: keep RDKit's coordinates
            for k, (x, y) in before.items():
                m.at[k]["x"], m.at[k]["y"] = x, y


def _arrange(frags, where, moves):
    """Order fragments so that partners of an arrow sit next to each other, turn each fragment (by 180°,
    never mirrored: a mirror image of a wedged drawing is the enantiomer) so its reacting atom faces its
    partner, and line the reacting atoms up vertically. Returns (ordered frags, x offsets, y offsets, width)."""
    by_name = {fr.name: fr for fr in frags}
    links = {}  # (name_a, name_b) -> list of (key_a, key_b)
    for mv in moves:
        atoms = [where[m] for m in mv.src + mv.dst if m in where]
        for i in range(len(atoms)):
            for j in range(i + 1, len(atoms)):
                (na, ka), (nb, kb) = atoms[i], atoms[j]
                if na != nb:
                    links.setdefault((na, nb), []).append((ka, kb))
                    links.setdefault((nb, na), []).append((kb, ka))
    partners = {fr.name: {b for (a, b) in links if a == fr.name} for fr in frags}
    radical = {
        fr.name
        for mv in moves
        if mv.electrons == 1
        for m in mv.src + mv.dst
        if m in where
        for fr in frags
        if fr.name == where[m][0]
    }
    # order: start from the fragment with most partners, grow a chain; spectators last
    order = []
    active = [fr.name for fr in frags if partners[fr.name]]
    if active:
        cur = max(active, key=lambda n: (len(partners[n]), len(by_name[n].m.at)))
        order.append(cur)
        while True:
            nxt = [n for n in active if n not in order and partners[n] & set(order)]
            if not nxt:
                rest = [n for n in active if n not in order]
                if not rest:
                    break
                nxt = rest
            n = nxt[0]
            # attach on the side of the partner already in the chain
            ends = (order[0], order[-1])
            if ends[0] in partners[n] and ends[1] not in partners[n]:
                order.insert(0, n)
            else:
                order.append(n)
    order += [fr.name for fr in frags if fr.name not in order]
    # face partners: a reacting atom with its partner to the right should be on the right half, and vice versa
    pos = {n: i for i, n in enumerate(order)}
    for n in order:
        m = by_name[n].m
        x0, _, x1, _ = _bbox(m)
        cx = (x0 + x1) / 2
        want = 0.0
        for (a, b), pairs in links.items():
            if a != n:
                continue
            d = 1 if pos[b] > pos[a] else -1
            for ka, _ in pairs:
                want += d * (m.at[ka]["x"] - cx)
        if want < -0.3:
            m.rot(180, cx, 0)
    # x offsets left to right; y: centre, then pull partners to the height of the reacting atom
    xs, ys, x = [], [], 0.0
    yoff = {}
    for n in order:
        m = by_name[n].m
        x0, y0, x1, y1 = _bbox(m)
        m.shift(-x0, -(y0 + y1) / 2)
    for i, n in enumerate(order):
        m = by_name[n].m
        dy = 0.0
        if i > 0 and (order[i - 1], n) in links:
            ka, kb = links[(order[i - 1], n)][0]
            prev = by_name[order[i - 1]].m
            dy = (prev.at[ka]["y"] + yoff.get(order[i - 1], 0.0)) - m.at[kb]["y"]
            dy = max(-1.5, min(1.5, dy))
        yoff[n] = dy
        if i > 0 and (order[i - 1], n) in links:
            # reacting atoms far enough apart for the arrows between them: two fish-hooks meet in the middle
            ka, kb = links[(order[i - 1], n)][0]
            prev_x = xs[i - 1] + by_name[order[i - 1]].m.at[ka]["x"]
            need = 2.8 if (n in radical and order[i - 1] in radical) else 2.2
            x = max(x, prev_x + need - m.at[kb]["x"])
        x0, _, x1, _ = _bbox(m)
        xs.append(x)
        ys.append(dy)
        x += (x1 - x0) + GAP_FRAG
    return [by_name[n] for n in order], xs, ys, x - GAP_FRAG


def render(mech, path=None):
    """Draw every state of the mechanism in a row. Returns the Canvas (its .issues() lists all violations)."""
    states = mech.states()
    panels = []  # (frags in order, where, width, x offsets, y offsets, y range)
    for p, smi in enumerate(states):
        frags, where = state_fragments(smi, f"s{p}")
        for fr in frags:
            _stretch_label_bonds(fr.m)
        moves = mech.steps[p].moves if p < len(mech.steps) else []
        frags, xs, ys, w = _arrange(frags, where, moves)
        ylo = min(_bbox(fr.m)[1] + dy for fr, dy in zip(frags, ys))
        yhi = max(_bbox(fr.m)[3] + dy for fr, dy in zip(frags, ys))
        panels.append((frags, where, w, xs, ys, ylo, yhi))
    width = sum(pn[2] for pn in panels) + GAP_PANEL * (len(panels) - 1) + 2 * PAD
    ylo = min(pn[5] for pn in panels)
    yhi = max(pn[6] for pn in panels)
    height = (yhi - ylo) + 2 * PAD + 1.5
    S = Canvas(int(width * B), int(height * B), B)
    cy = (PAD + 0.75 - ylo) * B  # pixel y of y = 0
    x = PAD
    placed = []
    for p, (frags, where, w, xs, ys, *_) in enumerate(panels):
        for fr, dx, dy in zip(frags, xs, ys):
            S.mol(fr.name, fr.m, (x + dx) * B, cy + dy * B)
        placed.append((frags, where))
        if p < len(panels) - 1:
            S.arrow((x + w + 0.6) * B, cy, (x + w + GAP_PANEL - 0.6) * B, cy)
        x += w + GAP_PANEL
    # unpaired electrons (facing the partner when the electron is used by an arrow) and charges on bare vertices
    for p, (frags, where) in enumerate(placed):
        moves = mech.steps[p].moves if p < len(mech.steps) else []
        aim = {}
        for mv in moves:
            if mv.electrons == 1 and len(mv.src) == 1 and mv.src[0] in where:
                other = [a for a in mv.dst if a != mv.src[0]]
                if other and other[0] in where:
                    aim[where[mv.src[0]]] = where[other[0]]
        for fr in frags:
            for k, a in fr.m.at.items():
                if a.get("_rad"):
                    prefer = None
                    if (fr.name, k) in aim:
                        tx, ty = S.P(*aim[(fr.name, k)])
                        px, py = S.P(fr.name, k)
                        prefer = math.degrees(math.atan2(-(ty - py), tx - px))
                    S.rad(fr.name, k, _free_angle(fr.m, k, prefer=prefer), gap=0.42 if a.get("label") else 0.36)
                if a.get("_q") and not a.get("label"):
                    ang = math.radians(_free_angle(fr.m, k, prefer=60))
                    px, py = S.P(fr.name, k)
                    S.chem(
                        px + math.cos(ang) * 0.42 * B, py - math.sin(ang) * 0.42 * B + 5, "+" if a["_q"] > 0 else "−"
                    )
    # curved arrows of each step on the state before it
    for k, st in enumerate(mech.steps):
        frags, where = placed[k]
        by_name = {fr.name: fr for fr in frags}
        for i, mv in enumerate(st.moves):
            src, dst, kw = _arrow_for(mv, where, by_name)
            hint = st.hints[i] if i < len(st.hints) else {}
            if src[0] == "lp" and (src[1], src[2]) not in S.__dict__.get("lps", {}):
                target = S.P(dst[1], dst[2])
                px, py = S.P(src[1], src[2])
                prefer = math.degrees(math.atan2(-(target[1] - py), target[0] - px))
                S.lone(src[1], src[2], _free_angle(by_name[src[1]].m, src[2], prefer=prefer))
            S.push(src, dst, side=hint.get("side"), **kw)
    if path:
        S.save(path)
    return S
