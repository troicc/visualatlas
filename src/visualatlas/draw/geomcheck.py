"""Generic geometry checks for drawn structural formulae.
Rules (IUPAC GR-1.1 / GR-1.10, ACS 1996 settings):
- Angle between neighbouring bonds at a non-metal atom: 120° in open chains (180° only for sp atoms,
  i.e. a triple bond or two double bonds); inside rings by ring size: 3 → 60°, 4 → 90°, 5 → 108°, 6 → 120°;
  an exocyclic bond on a ring atom bisects the exterior angle.
- Four-coordinate (wedged) centres: every neighbouring angle within 50°–150°; with one bold and one hashed
  wedge the IUPAC ST-1.1.1 template (60° / 90° / 120°) is enforced.
- Bond length: carbon skeleton 1.0, both ends labelled 1.0–1.5, metal–ligand 1.0–2.2.
Atoms marked persp (bridged rings in perspective, 3D views) are reported but not judged."""

import math
import re
from collections import deque

METALS = {
    "Rh",
    "Ru",
    "Ti",
    "Os",
    "Pd",
    "Pt",
    "Ir",
    "Co",
    "Ni",
    "Fe",
    "Cu",
    "Zn",
    "Mn",
    "Cr",
    "Mo",
    "W",
    "Al",
    "Li",
    "Na",
    "K",
    "Mg",
    "Cs",
    "Ag",
    "Au",
    "Sc",
    "Zr",
    "Hf",
    "V",
}
SKIP = ("d", "dot", "arc", "hp", "free")  # 'p'（部分双键）算几何
TET = "ST-1.1.1"


def is_metal(lab):
    if not lab:
        return False
    m = re.match(r"\[?([A-Z][a-z]?)", lab)
    return bool(m) and m.group(1) in METALS


def ring_path(adj, k, u, v, maxn=8):
    """不经过 k 的 u→v 最短路（含两端）"""
    if u == v:
        return None
    prev = {u: None}
    q = deque([u])
    while q:
        x = q.popleft()
        for y in adj[x]:
            if y == k or y in prev:
                continue
            prev[y] = x
            if y == v:
                path = [v]
                while prev[path[-1]] is not None:
                    path.append(prev[path[-1]])
                return path if len(path) + 1 <= maxn else None
            q.append(y)
    return None


def ring_size(adj, k, u, v, maxn=8):
    """u、v 都是 k 的邻居：不经过 k 的 u→v 最短路 + k 构成的环大小"""
    if u == v:
        return None
    seen = {u: 0}
    q = deque([u])
    while q:
        x = q.popleft()
        if seen[x] >= maxn:
            break
        for y in adj[x]:
            if y == k or y in seen:
                continue
            seen[y] = seen[x] + 1
            if y == v:
                return seen[y] + 2
            q.append(y)
    return None


def check_generic(name, atoms, bonds, persp=(), tsring=()):
    """atoms: {k: (x, y, label)}（y 朝上或朝下都行）；bonds: [(a, b, t)]"""
    out = []
    adj = {k: [] for k in atoms}
    bl = [(a, b, t) for a, b, t in bonds if t not in SKIP or t == "d"]  # 正在形成的虚线键也算几何（过渡态里的四面体碳）
    for a, b, t in bl:
        adj[a].append(b)
        adj[b].append(a)
    order = {}
    for a, b, t in bl:
        order[(a, b)] = order[(b, a)] = t
    for k, nbs in adj.items():
        x0, y0, lab = atoms[k]
        if is_metal(lab) or len(nbs) < 2:
            continue
        angs = sorted(((math.degrees(math.atan2(atoms[o][1] - y0, atoms[o][0] - x0)) % 360), o) for o in nbs)
        n = len(angs)
        gaps = [((angs[(i + 1) % n][0] - angs[i][0]) % 360, angs[i][1], angs[(i + 1) % n][1]) for i in range(n)]
        if n == 2:
            g = min(gaps[0][0], 360 - gaps[0][0])
            gaps = [(g, gaps[0][1], gaps[0][2])]
        ts = [order[(k, o)] for o in nbs]
        linear = 3 in ts or ts.count(2) == 2
        tag = f"{name}: {k}({lab or 'C'})"
        # 三配位以上的原子：任何两根键（不一定相邻）都不许共线——共线就是“炔烃式”直线（2001 期 Ru 图里的 C–H 楔与 C–CH2）
        if n >= 3 and k not in persp:
            for i_ in range(n):
                for j_ in range(i_ + 1, n):
                    g2 = abs(angs[i_][0] - angs[j_][0]) % 360
                    g2 = min(g2, 360 - g2)
                    if g2 >= 168:
                        out.append(f"{tag} {angs[i_][1]}–{k}–{angs[j_][1]} 共线 {g2:.0f}°")
        if (
            n == 4
            and sorted(ts, key=str).count("w") == 1
            and ts.count("h") == 1
            and not any(ring_size(adj, k, p, q) for i, p in enumerate(nbs) for q in nbs[i + 1 :])
        ):
            # IUPAC ST-1.1.1：实楔与虚楔相邻 60°，两根平面键 120°，其余 90°
            kind = {o: order[(k, o)] for o in nbs}
            bad = []
            for g, u, v in gaps:
                pair = {kind[u], kind[v]}
                if pair == {"w", "h"}:
                    want = 60
                elif "w" not in pair and "h" not in pair:
                    want = 120
                else:
                    want = 90
                if abs(g - want) > 12:
                    bad.append(f"{u}–{k}–{v} = {g:.0f}°（模板 {want}°）")
            if bad:
                out.append(f"{tag} 四面体模板不符：" + "；".join(bad))
            continue
        for g, u, v in gaps:
            rs = ring_size(adj, k, u, v)
            if n == 2 and linear:
                want, tol = 180, 6
            elif rs and rs <= 6:
                want, tol = 180 - 360 / rs, 8
                if order[(k, u)] == "d" or order[(k, v)] == "d":
                    tol = 16  # 含正在形成的键的过渡态环：放宽
                pth = ring_path(adj, k, u, v)
                if pth and any(is_metal(atoms[x][2]) for x in pth):
                    tol = 16  # 螯合环：M–L 键长 1.5，环不可能正多边形
            elif rs:
                want, tol = None, None  # 大环：只查不太离谱
            elif n <= 3:
                # 环外键：若 k 在环上，平分外角
                ring_int = [
                    180 - 360 / s
                    for s in (ring_size(adj, k, p, q) for i, p in enumerate(nbs) for q in nbs[i + 1 :])
                    if s and s <= 6
                ]
                want = (360 - ring_int[0]) / 2 if (ring_int and n == 3) else 120
                tol = 8
            else:
                want, tol = None, None
            if k in tsring and n == 2:
                if not 95 <= g <= 135:
                    out.append(f"{tag} {u}–{k}–{v} = {g:.0f}°（过渡态环 95–135°）")
                continue
            if want is None:
                if g < 50 or g > 150 and n >= 4 or g > 190:
                    out.append(f"{tag} {u}–{k}–{v} = {g:.0f}°")
                continue
            if abs(g - want) > tol:
                out.append(f"{tag} {u}–{k}–{v} = {g:.0f}°（应 {want:.0f}°{'，透视' if k in persp else ''}）")
    for a, b, t in bl:
        if t == "d":
            continue
        L = math.dist(atoms[a][:2], atoms[b][:2])
        la, lb = atoms[a][2], atoms[b][2]
        if is_metal(la) or is_metal(lb):
            ok = 0.95 <= L <= 2.2
        elif la and lb:
            ok = 0.96 <= L <= 1.54
        else:
            ok = abs(L - 1) <= 0.04
        if not ok:
            out.append(f"{name}: 键 {a}–{b} 长 {L:.2f}")
    return out


def collinear(name, atoms, bonds, persp=(), min_n=3):
    """三配位以上的非金属原子：任何两根键共线（≥ 168°）都报——sp 原子除外"""
    out = []
    adj = {k: [] for k in atoms}
    order = {}
    for a, b, t in bonds:
        if t in SKIP and t != "d":
            continue
        adj[a].append(b)
        adj[b].append(a)
        order[(a, b)] = order[(b, a)] = t
    for k, nbs in adj.items():
        x0, y0, lab = atoms[k]
        if is_metal(lab) or lab == "M" or len(nbs) < min_n or k in persp:
            continue
        if min_n == 2 and len(nbs) == 2 and not any(is_metal(atoms[o][2]) or atoms[o][2] == "M" for o in nbs):
            continue  # 纯有机二配位已由 120° 规则管
        ts = [order[(k, o)] for o in nbs]
        if 3 in ts or ts.count(2) >= 2:
            continue
        if lab == "H" and len(nbs) == 2 and "d" in ts:
            continue  # 正在转移的 H（X–H···Y）本来就接近直线
        angs = [(math.degrees(math.atan2(atoms[o][1] - y0, atoms[o][0] - x0)) % 360, o) for o in nbs]
        for i in range(len(angs)):
            for j in range(i + 1, len(angs)):
                g = abs(angs[i][0] - angs[j][0]) % 360
                g = min(g, 360 - g)
                if g >= 168:
                    out.append(f"{name}: {k}({lab or 'C'}) {angs[i][1]}–{k}–{angs[j][1]} 共线 {g:.0f}°")
    return out
