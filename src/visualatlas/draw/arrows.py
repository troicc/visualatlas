"""Curved (electron-pushing) arrows: placement by exhaustive search under hard geometric checks.

Every arrow is a true circular arc. Its shape is the arch ratio k = arch height / chord:
0.18–0.42 for ordinary arrows, 0.40–0.70 for "own" curls (bond electrons collapsing onto one of its
own atoms, or a lone pair forming a bond on its own atom; k > 0.5 is more than a semicircle).
Symmetric cubic Béziers cannot draw the C-shaped curl (C=O → O, N: → N–C) and are not used.

Hard checks on every candidate (a candidate that fails any of them is rejected):
  - the tangent at the tip points at the target (atom ≤ 25°, bond ≤ 50°);
  - the tip is closer to its target than to any other atom or glyph, and does not enter along a bond of the target;
  - the shaft crosses no bond, passes through no label, and keeps clear of bonds, labels, lone-pair dots,
    text and other arrows;
  - the tail leaves a lone pair outwards, and starts 0.22 B off a bond centre line (+0.12 B for a double bond);
  - chord length 0.75–3.2 B (curls may be shorter);
  - side='up'|'down'|'left'|'right' forces the bulge direction.
If no candidate passes, the best-scoring one is drawn and an issue is recorded on the canvas.

Endpoints:
  src: ('lp', mol, atom) | ('bond', mol, a, b[, u]) | ('rad', mol, atom)
  dst: ('atom', mol, atom) | ('bond', mol, a, b) | ('form', mol, a, b) | ('mid', molA, a, molB, b)
Usage: push(canvas, src, dst, side=None, tip_dir=None, u=None, k=None, half=False)
  half=True draws a fish-hook (single-electron) arrow.
"""

import math

import numpy as np

from .mol import K, label_box, layout, adv

ACC = "#B3261E"


def f(v):
    return f"{v:.1f}"


SIDES = {"up": (0, -1), "down": (0, 1), "left": (-1, 0), "right": (1, 0)}


# ---------------- 几何基础 ----------------
def _arc_pts(p0, p1, k, s, n=64):
    """从 p0 到 p1 的圆弧，拱高 = k·弦长，s = ±1 选拱向（屏幕坐标下弦左侧为 +）。返回 (n, 2) 点列"""
    p0 = np.array(p0, float)
    p1 = np.array(p1, float)
    d = p1 - p0
    c = np.hypot(*d)
    if c < 1e-6:
        return np.repeat(p0[None], n, 0)
    nrm = np.array([d[1], -d[0]]) / c * s  # 拱出去的方向
    t = k * c
    h = c / 2
    R = (h * h + t * t) / (2 * t)
    mid = (p0 + p1) / 2
    ctr = mid + nrm * (t - R)
    a0 = math.atan2(p0[1] - ctr[1], p0[0] - ctr[0])
    a1 = math.atan2(p1[1] - ctr[1], p1[0] - ctr[0])
    bul = mid + nrm * t
    am = math.atan2(bul[1] - ctr[1], bul[0] - ctr[0])

    def norm(a):
        return (a + 2 * math.pi) % (2 * math.pi)

    # 选择经过拱顶的那一圈方向
    ccw = norm(am - a0) < norm(a1 - a0)
    sweep = norm(a1 - a0) if ccw else -norm(a0 - a1)
    ang = a0 + sweep * np.linspace(0, 1, n)
    return np.stack([ctr[0] + R * np.cos(ang), ctr[1] + R * np.sin(ang)], 1)


def _seg_dist(P, A, Bq):
    """点集 P (n,2) 到线段集 (m,2),(m,2) 的距离矩阵 (n,m)"""
    if len(A) == 0:
        return np.full((len(P), 0), 1e9)
    v = Bq - A
    L2 = (v**2).sum(1)
    L2[L2 == 0] = 1
    w = P[:, None, :] - A[None, :, :]
    t = np.clip((w * v[None]).sum(2) / L2[None], 0, 1)
    q = A[None] + t[..., None] * v[None]
    return np.hypot(*(P[:, None, :] - q).transpose(2, 0, 1))


def _box_dist(P, boxes):
    """点集到矩形集的距离（矩形内为负的穿入深度近似为 0 以下）"""
    if len(boxes) == 0:
        return np.full((len(P), 0), 1e9)
    b = np.array(boxes, float)
    dx = np.maximum(np.maximum(b[None, :, 0] - P[:, None, 0], 0), P[:, None, 0] - b[None, :, 2])
    dy = np.maximum(np.maximum(b[None, :, 1] - P[:, None, 1], 0), P[:, None, 1] - b[None, :, 3])
    inside = (dx == 0) & (dy == 0)
    d = np.hypot(dx, dy)
    d[inside] = -1
    return d


def _cross(P, A, Bq):
    """折线 P 是否与线段集相交：返回布尔 (m,)"""
    if len(A) == 0:
        return np.zeros(0, bool)
    a = P[:-1]
    b = P[1:]

    def orient(p, q, r):
        return (q[..., 0] - p[..., 0]) * (r[..., 1] - p[..., 1]) - (q[..., 1] - p[..., 1]) * (r[..., 0] - p[..., 0])

    a_ = a[:, None, :]
    b_ = b[:, None, :]
    c_ = A[None]
    d_ = Bq[None]
    d1 = orient(a_, b_, c_)
    d2 = orient(a_, b_, d_)
    d3 = orient(c_, d_, a_)
    d4 = orient(c_, d_, b_)
    return ((d1 * d2 < 0) & (d3 * d4 < 0)).any(0)


# ---------------- 场景：所有障碍物 ----------------
def _scene(S):
    B = S.B
    bonds = []
    labels = []
    for name, (m, ox, oy) in S.mols.items():
        for e in m.bd:
            if e.get("t") == "arc":
                continue
            pa = (ox + m.at[e["a"]]["x"] * B, oy + m.at[e["a"]]["y"] * B)
            pb = (ox + m.at[e["b"]]["x"] * B, oy + m.at[e["b"]]["y"] * B)
            bonds.append((name, e["a"], e["b"], e.get("t"), pa, pb))
        for k, a in m.at.items():
            if a.get("label"):
                bx = label_box(a, B)
                labels.append((name, k, (bx[0] + ox, bx[1] + oy, bx[2] + ox, bx[3] + oy)))
    texts = [bx for _, bx in S.__dict__.get("tboxes", [])]
    chg = S.__dict__.get("chgboxes", [])
    dots = [
        (n, k, (c[0] - 0.09 * B, c[1] - 0.14 * B, c[0] + 0.09 * B, c[1] + 0.14 * B))
        for (n, k), (c, _) in S.__dict__.get("lps", {}).items()
    ]
    curls = [c["pts"] for c in S.__dict__.get("curls", [])]
    paths = [np.array(p, float) for _, p in S.__dict__.get("apaths", []) if len(p) > 1]
    return bonds, labels, texts, chg, dots, curls, paths


def _atom_xy(S, name, k):
    m, ox, oy = S.mols[name]
    a = m.at[k]
    return np.array([ox + a["x"] * S.B, oy + a["y"] * S.B])


def _label_box(S, name, k):
    m, ox, oy = S.mols[name]
    a = m.at[k]
    if not a.get("label"):
        return None
    b = label_box(a, S.B)
    return (b[0] + ox, b[1] + oy, b[2] + ox, b[3] + oy)


def _glyph_box(S, name, k):
    """锚定字符本身的框（多字形缩写如 ²⁻O₃P 时，箭头尖要对准 P 这个字，而不是整串标签的边，第二轮 PGM 评审）"""
    m, ox, oy = S.mols[name]
    a = m.at[k]
    if not a.get("label"):
        return None
    B = S.B
    fs = K["fs"] * B
    mg = K["mg"] * B
    g = [q for q in layout(a["label"], fs, a.get("anchor", 0)) if q[5] is not None and q[2] == 0]
    i = min(a.get("anchor", 0), len(g) - 1) if g else 0
    if not g:
        return _label_box(S, name, k)
    q = g[i]
    x = ox + a["x"] * B
    y = oy + a["y"] * B
    return (x + q[1] - mg, y - fs * 0.37 - mg, x + q[1] + adv(q[0]) * q[3] + mg, y + fs * 0.37 + mg)


def _ray_exit(c, box, u):
    """从 c 沿单位向量 u 出发离开矩形 box 的距离"""
    ts = []
    for i, (lo, hi) in enumerate(((box[0], box[2]), (box[1], box[3]))):
        if abs(u[i]) > 1e-9:
            ts.append(((hi if u[i] > 0 else lo) - c[i]) / u[i])
    return max(0, min(ts)) if ts else 0


def _bond_of(S, name, a, b):
    m = S.mols[name][0]
    for e in m.bd:
        if {e["a"], e["b"]} == {a, b}:
            return e
    return None


# ---------------- 箭尾与箭头尖的候选 ----------------
def _tails(S, src):
    B = S.B
    if src[0] == "lp":
        lps = S.__dict__.get("lps", {})
        if (src[1], src[2]) in lps:
            c, u = lps[(src[1], src[2])]
            c = np.array(c)
            u = np.array(u)
            return [dict(p=c + u * 0.22 * B, dir=u, kind="lp")]  # 箭尾离两点有肉眼可见的空隙（第二轮：0.15 B 看着贴点）
        # 没画孤对点：从原子四周 8 个方向出发（不推荐，L8 式图应先 S.lone）
        xy = _atom_xy(S, src[1], src[2])
        box = _label_box(S, src[1], src[2])
        out = []
        for th in range(0, 360, 45):
            u = np.array([math.cos(math.radians(th)), -math.sin(math.radians(th))])
            d = (_ray_exit(xy, box, u) if box else 0.2 * B) + 0.15 * B
            out.append(dict(p=xy + u * d, dir=u, kind="lp"))
        return out
    if src[0] == "rad":
        p = np.array(S.rads[(src[1], src[2])])
        return [dict(p=p, dir=None, kind="rad")]
    # 键：中点附近、线外 0.22 B（双键 +0.12 B），两侧都试
    A = _atom_xy(S, src[1], src[2])
    Bq = _atom_xy(S, src[1], src[3])
    v = Bq - A
    L = np.hypot(*v)
    n = np.array([-v[1], v[0]]) / L
    e = _bond_of(S, src[1], src[2], src[3])
    dbl = e is not None and e["t"] in (2, "p", 3)
    us = (src[4],) if len(src) > 4 and src[4] is not None else (0.4, 0.5, 0.6)
    out = []
    for u in us:
        for sg in (1, -1):
            off = (0.22 + (0.12 if dbl else 0)) * B
            out.append(dict(p=A + v * u + n * sg * off, dir=n * sg, kind="bond", u=u))
    return out


def _tips(S, dst, tail, tip_dir=None, tip_win=32):
    B = S.B
    if dst[0] == "atom":
        c = _atom_xy(S, dst[1], dst[2])
        box = _glyph_box(S, dst[1], dst[2])
        out = []
        if box:
            c = np.array([(box[0] + box[2]) / 2, (box[1] + box[3]) / 2])
        rng = range(0, 360, 12)
        for th in rng:
            if tip_dir is not None and min(abs(th - tip_dir) % 360, 360 - abs(th - tip_dir) % 360) > tip_win:
                continue
            u = np.array([math.cos(math.radians(th)), -math.sin(math.radians(th))])
            mg = K["mg"] * B
            gb = (box[0] + mg, box[1] + mg, box[2] - mg, box[3] - mg) if box else None
            d = (
                (_ray_exit(c, gb, u) + 0.12 * B) if box else 0.28 * B
            )  # 第三轮：尖到原子中心 0.38 B 左右（原 0.65 B 偏远）
            out.append(dict(p=c + u * d, target=c, kind="atom", th=th))
        return out
    if dst[0] in ("bond", "form"):
        A = _atom_xy(S, dst[1], dst[2])
        Bq = _atom_xy(S, dst[1], dst[3])
        v = Bq - A
        L = np.hypot(*v)
        n = np.array([-v[1], v[0]]) / L
        e = _bond_of(S, dst[1], dst[2], dst[3])
        dbl = e is not None and e["t"] in (2, "p", 3)
        same = 1 if np.dot(tail["p"] - (A + v * 0.5), n) >= 0 else -1
        out = []
        # 两侧都试（环上 π 推移时，箭尾在环外而按整条直线算「同侧」会落进环内，第二轮 GAPDH）；
        # 横穿目标键由相交检查拦住，同侧的另加分
        for sg in (same, -same):
            for u in (0.42, 0.5, 0.58):
                off = (0.28 + (0.1 if dbl else 0)) * B
                out.append(dict(p=A + v * u + n * sg * off, target=A + v * u, kind="bond", n=n * sg, same=(sg == same)))
        return out
    if dst[0] == "mid":
        A = _atom_xy(S, dst[1], dst[2])
        Bq = _atom_xy(S, dst[3], dst[4])
        m = (A + Bq) / 2
        return [dict(p=m, target=m, kind="mid")]
    raise ValueError(dst)


# ---------------- 评估一个候选 ----------------
def _evaluate(S, src, dst, tail, tip, k, s, scene, side=None, own=False, hint=0):
    B = S.B
    bonds, labels, texts, chg, dots, curls, paths = scene
    p0, p1 = tail["p"], tip["p"]
    chord = np.hypot(*(p1 - p0))
    bad = []
    lo = (
        (0.26 if (dst[0] == "atom" and _glyph_box(S, dst[1], dst[2]) is None) else 0.48) if own else 0.75
    )  # 第四轮第二组：卷弧下限 0.36 → 0.48 B（评审：0.45 B 的卷弧肉眼就是小钩）   # 第四轮：键 → 自己一端隐式顶点的卷弧放宽到 0.26 B
    if (
        own
        and src[0] == "bond"
        and dst[0] == "atom"
        and _glyph_box(S, dst[1], dst[2]) is None
        and any(_label_box(S, src[1], kk) for kk in src[2:4] if kk != dst[2])
    ):
        lo = 0.30  # 另一端是单字母标签（C–H → C）：弦长可低到 0.30 B，但可见箭身 ≥ 0.45 B、拱高 ≥ 0.35 弦长（r4/rad 第三轮评审：0.20 B 放过了缩成小点的鱼钩）      # 卷弧（1.0 B 短键、环内孤对上的 L8 小弧）下限 0.36 B；普通箭头 0.75 B
    if chord < lo * B:
        bad.append(f"弦长 {chord / B:.2f} B 太短（{'卷弧' if own else '普通箭头'}下限 {lo} B）")
    if chord > 3.2 * B:
        bad.append(f"弦长 {chord / B:.2f} B 太长（> 3.2 B），应把反应伙伴挪近")
    P = _arc_pts(p0, p1, k, s)
    n = len(P)
    uu = np.linspace(0, 1, n)
    # 拱向硬约束
    mid = (p0 + p1) / 2
    bul = P[n // 2] - mid
    if side:
        sv = np.array(SIDES[side])
        nb = np.hypot(*bul) or 1
        if np.dot(bul / nb, sv) < 0.35:
            bad.append(f"拱向不是 {side}")
    # 尖的方向：末端切线指向目标
    tang = _head_dir(P, 0.27 * B)
    to = tip["target"] - p1
    dto = np.hypot(*to)
    if tip["kind"] in ("atom", "bond") and dto > 1e-6:
        cosang = float(np.dot(tang, to / dto))
        lim = math.cos(math.radians(25 if tip["kind"] == "atom" else 50))
        if cosang < lim:
            bad.append(f"箭头没有朝向目标（偏 {math.degrees(math.acos(max(-1, min(1, cosang)))):.0f}°）")
    else:
        cosang = 1.0
    # 箭尾方向：从孤对出发时先离开原子
    if tail["kind"] == "lp" and tail["dir"] is not None:
        t0 = P[3] - P[0]
        t0 = t0 / (np.hypot(*t0) or 1)
        if np.dot(t0, tail["dir"]) < math.cos(math.radians(80)):
            bad.append("箭尾没有从孤对点向外出发")
    # 障碍物
    A = np.array([b[4] for b in bonds]) if bonds else np.zeros((0, 2))
    Bq = np.array([b[5] for b in bonds]) if bonds else np.zeros((0, 2))
    crossed = _cross(P, A, Bq)
    if crossed.any():
        i = int(np.argmax(crossed))
        bad.append(f"箭身与键 {bonds[i][0]}.{bonds[i][1]}–{bonds[i][2]} 相交")
    D = _seg_dist(P, A, Bq)
    src_bond = (src[1], src[2], src[3]) if src[0] == "bond" else None
    clear = 9e9
    for j, b in enumerate(bonds):
        key = {b[1], b[2]}
        isrc = src_bond and b[0] == src_bond[0] and key == {src_bond[1], src_bond[2]}
        itgt = dst[0] in ("bond", "form") and b[0] == dst[1] and key == {dst[2], dst[3]}
        near_src_atom = (src[0] in ("lp", "rad") and b[0] == src[1] and src[2] in key) or (
            src_bond and b[0] == src_bond[0] and (key & {src_bond[1], src_bond[2]})
        )
        near_tgt_atom = dst[0] == "atom" and b[0] == dst[1] and dst[2] in key
        lim = 0.14 * B + (0.1 * B if b[3] in (2, "p") else 0)
        lo = 0.3 if near_src_atom else 0.15
        hi = 0.7 if near_tgt_atom else 0.85
        sel = (uu > lo) & (uu < hi)
        if not (isrc or itgt) and sel.any() and D[sel, j].min() < lim:
            bad.append(f"箭身贴着键 {b[0]}.{b[1]}–{b[2]}（{D[sel, j].min() / B:.2f} B）")
        # 两端附近也不能真的压在键上（≥ 0.06 B）
        if not (isrc or itgt) and D[:, j].min() < 0.06 * B:
            bad.append(f"箭身压在键 {b[0]}.{b[1]}–{b[2]} 上")
        # 箭头三角区离非目标键 ≥ 0.2 B
        if not itgt and not near_tgt_atom:
            dh = D[uu > 0.9, j].min()
            if dh < 0.2 * B + (0.08 * B if b[3] in (2, "p") else 0):
                bad.append(f"箭头尖贴着键 {b[0]}.{b[1]}–{b[2]}（{dh / B:.2f} B）")
        if not (isrc or itgt):
            clear = min(clear, D[(uu > 0.2) & (uu < 0.8), j].min())
    # 标签
    LB = [l[2] for l in labels]
    Dl = _box_dist(P, LB)
    for j, l in enumerate(labels):
        own_src = (src[0] in ("lp", "rad") and (l[0], l[1]) == (src[1], src[2])) or (
            src_bond and l[0] == src_bond[0] and l[1] in (src_bond[1], src_bond[2])
        )
        own_tgt = (dst[0] == "atom" and (l[0], l[1]) == (dst[1], dst[2])) or (
            dst[0] in ("bond", "form") and l[0] == dst[1] and l[1] in (dst[2], dst[3])
        )
        dl = Dl[:, j]
        pen = dl < 0
        if own_tgt:
            pen &= uu < 0.9  # 尖本身停在字形外 0.12 B，可以进入目标标签的留白区（第三轮）
        if pen.any():
            bad.append(f"箭身穿过标签 {l[0]}.{l[1]}")
            continue
        sel = np.ones(n, bool)
        if own_src:
            sel &= uu > 0.25
        if own_tgt:
            sel &= uu < 0.75
        if sel.any() and dl[sel].min() < 0.12 * B:
            bad.append(f"箭身贴着标签 {l[0]}.{l[1]}（{dl[sel].min() / B:.2f} B）")
        if sel.any():
            clear = min(clear, dl[sel].min())
    for bx in texts + chg:
        dl = _box_dist(P, [bx])[:, 0]
        if (dl[(uu > 0.05) & (uu < 0.97)] < 0.1 * B).any():
            bad.append("箭身贴着 / 穿过文字或电荷")
            break
        clear = min(clear, dl.min())
    for nm, kk, bx in dots:
        if src[0] == "lp" and (nm, kk) == (src[1], src[2]):
            continue
        dl = _box_dist(P, [bx])[:, 0]
        if dl.min() < 0.1 * B:
            bad.append(f"箭身贴着 {nm}.{kk} 的孤对点")
            break
    for j, Q in enumerate(curls + paths):
        Pq = P
        if dst[0] == "mid" and j < len(curls):
            # two fish-hooks that make one bond meet between its atoms by design: compare only the first 75 %
            # of both shafts, so the tips may come close (added 2026-11-07, not in the ported engine)
            other = S.__dict__["curls"][j]["dst"]
            if other[0] == "mid" and {tuple(other[1:3]), tuple(other[3:5])} == {tuple(dst[1:3]), tuple(dst[3:5])}:
                Pq, Q = P[: int(len(P) * 0.75)], Q[: int(len(Q) * 0.75)]
        dq = np.hypot(*(Pq[:, None, :] - Q[None]).transpose(2, 0, 1)).min()
        if dq < 0.22 * B:
            bad.append("与另一根箭头相交或贴近")
            break
        clear = min(clear, dq)
    # 尖不能顺着目标原子的某根键进来（第二轮：C6=N1→N 的尖离 R–N 键比离 N 还近；H₃C–C 键贴着负氢箭头的三角）
    if tip["kind"] == "atom":
        tv = p1 - tip["target"]
        tvn = np.hypot(*tv) or 1
        for b in bonds:
            if b[0] != dst[1] or dst[2] not in (b[1], b[2]):
                continue
            other = np.array(b[5] if b[1] == dst[2] else b[4])
            bv = other - tip["target"]
            bn = np.hypot(*bv) or 1
            if np.dot(tv / tvn, bv / bn) > math.cos(math.radians(38)):
                bad.append(f"箭头尖顺着键 {b[1]}–{b[2]} 进入目标")
                break
    # 尖离目标字形比离任何别的字形近（多字形缩写 ²⁻O₃P：尖不能落在 O₃ 上方，第二轮 PGM）
    if tip["kind"] == "atom":
        dt = np.hypot(*(p1 - tip["target"]))
        for name, (m, ox, oy) in S.mols.items():
            for kk, a in m.at.items():
                if not a.get("label"):
                    continue
                fs = K["fs"] * B
                for q in layout(a["label"], fs, a.get("anchor", 0)):
                    if q[5] is None or q[2] != 0:
                        continue
                    gc = np.array([ox + a["x"] * B + q[1] + adv(q[0]) * q[3] / 2, oy + a["y"] * B])
                    if np.hypot(*(gc - tip["target"])) < 2:
                        continue
                    if np.hypot(*(p1 - gc)) < dt - 0.05 * B:
                        bad.append(f"箭头尖离 {name}.{kk} 的「{q[0]}」比离目标字形还近")
                        break
    # 尖的归属：离目标比离别的原子近
    if tip["kind"] == "atom":
        dt = np.hypot(*(p1 - tip["target"]))
        for name, (m, ox, oy) in S.mols.items():
            for kk, a in m.at.items():
                if (name, kk) == (dst[1], dst[2]):
                    continue
                q = np.array([ox + a["x"] * B, oy + a["y"] * B])
                if np.hypot(*(p1 - q)) < dt - 0.05 * B:
                    bad.append(f"箭头尖离 {name}.{kk} 比离目标还近")
                    break
    # 打分
    k0 = 0.55 if own else 0.3
    score = (
        min(clear, (0.45 if own else 0.7) * B)
        + 0.8 * B * cosang
        + (0.1 * B if tip.get("same", True) else 0)
        - 0.8 * B * abs(k - k0)
        + (0.05 * B if hint and s == hint else 0)
    )
    if tail["kind"] == "bond":
        score -= 0.3 * B * abs(tail["u"] - 0.5)
    return bad, score, P


def _head_dir(P, hl):
    """箭头方向 = 箭身裁掉箭头长度后的末端 → 尖（第二轮：紧卷弧上精确切线与看得见的箭身方向差 30–50°，箭头像折了）"""
    seg = np.hypot(*np.diff(P, axis=0).T)
    cum = np.concatenate([[0], np.cumsum(seg)])
    i = int(np.searchsorted(cum, cum[-1] - hl * 0.85))
    i = max(0, min(i, len(P) - 2))
    d = P[-1] - P[i]
    return d / (np.hypot(*d) or 1)


def _is_own(src, dst):
    """卷弧类：键电子落回自己一端；孤对 → 本原子的键；键 → 相邻的键"""
    if src[0] == "bond" and dst[0] == "atom" and src[1] == dst[1] and dst[2] in src[2:4]:
        return True
    if src[0] == "lp" and dst[0] in ("bond", "form") and src[1] == dst[1] and src[2] in dst[2:4]:
        return True
    if src[0] == "bond" and dst[0] in ("bond", "form") and src[1] == dst[1] and set(src[2:4]) & set(dst[2:4]):
        return True
    return False


def _draw(S, P, col, half=False):
    B = S.B
    hl = 0.27 * B
    # 沿弧长裁掉箭头长度
    seg = np.hypot(*np.diff(P, axis=0).T)
    cum = np.concatenate([[0], np.cumsum(seg)])
    Ltot = cum[-1]
    keep = cum <= Ltot - hl * 0.85
    Q = P[keep]
    if len(Q) < 2:
        Q = P[:2]
    end = P[-1]
    tan = P[-1] - Q[-1]
    ang = math.atan2(tan[1], tan[0])
    w = K["lw"] * B * 1.1
    d = f"M{f(Q[0][0])} {f(Q[0][1])} " + " ".join(f"L{f(x)} {f(y)}" for x, y in Q[1:])
    S.parts.append(
        f'<path d="{d}" fill="none" stroke="{col}" stroke-width="{f(w)}" stroke-linecap="round" stroke-linejoin="round"/>'
    )
    S.parts.append(S._head(end[0], end[1], ang, col, 0.9, half=(1 if half else 0)))


def push(
    S,
    src,
    dst,
    bend=None,
    col=ACC,
    f=None,
    gap=None,
    rot=None,
    auto=None,
    tip_dir=None,
    half=False,
    best=None,
    side=None,
    u=None,
    k=None,
    tip_win=32,
    **kw,
):
    B = S.B
    if u is not None and src[0] == "bond":
        src = tuple(src[:4]) + (u,)
    elif src[0] == "bond" and len(src) > 4 and src[4] in ("bulge",):
        src = tuple(src[:4])
    elif src[0] == "bond" and len(src) > 5:
        src = tuple(src[:5])
    own = _is_own(src, dst)
    ks = (k,) if k else ((0.4, 0.48, 0.56, 0.64, 0.72) if own else (0.18, 0.24, 0.3, 0.36, 0.42))
    hint = 0 if bend is None or bend == 0 else (1 if bend > 0 else -1)
    scene = _scene(S)
    best_ok = None
    best_bad = None
    for tl in _tails(S, src):
        for tp in _tips(S, dst, tl, tip_dir, tip_win):
            for kk in ks:
                for s in (1, -1):
                    bad, sc, P = _evaluate(S, src, dst, tl, tp, kk, s, scene, side, own, hint)
                    if not bad:
                        if best_ok is None or sc > best_ok[0]:
                            best_ok = (sc, P, kk, s)
                    else:
                        pen = sc - 2 * B * len(bad) - (5 * B if any(x.startswith("拱向") for x in bad) else 0)
                        if best_bad is None or pen > best_bad[0]:
                            best_bad = (pen, P, bad)
    if best_ok:
        P = best_ok[1]
    else:
        P = best_bad[1]
        S.__dict__.setdefault("arrow_issues", []).append(
            f"弯箭头 {src[1:4]}→{dst[1:4]}：无可行弧，最接近的方案仍有 {'；'.join(best_bad[2][:3])}——改结构摆位或给 side/tip_dir"
        )
    _draw(S, P, col, half)
    S.__dict__.setdefault("curls", []).append(dict(src=src, dst=dst, pts=P))
    # 兼容旧检查（bio._crowd_issues 读 S.arrows 的 (bend, p0, p1)）：登记为「已由 l8arrow 检查」
    S.__dict__.setdefault("arrows_l8", []).append((src, dst, P))
    return S
