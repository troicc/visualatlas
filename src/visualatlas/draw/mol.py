"""Drawing model for structural formulae: atoms with 2D coordinates and labels, bonds with a drawing type.

This is a *drawing* model, not a chemical one: an atom is a position plus an optional label string
("CH2OH", "Rh"), and a bond is a pair of atom keys plus a drawing type. Chemistry (elements, charges,
lone pairs) lives in the RDKit graph; this module only lays out and renders it.

Units: one bond length = 1, screen y points down, angles are mathematical (0 = right, 90 = up).
Line weights follow the ACS 1996 document settings: bond length 0.508 cm, line width 0.04 cm,
bold wedge 0.071 cm (base width taken as 0.2), double-bond spacing 18 %, label margin 0.056 cm,
hash spacing 0.088 cm, Arial 10 pt.

Stereo and coordination geometry follow IUPAC 2006 (Brecher, Pure Appl. Chem. 78, 1897):
  ST-0.3  the narrow end of a wedge is at the stereocentre; dashed lines are only for partial / H bonds;
  ST-2.6  square planar: two bold wedges below (angle > 90°), two hashed wedges above, mirror-symmetric;
  ST-2.11 octahedral: two straight axial bonds, two bold wedges below, two hashed wedges above;
  GR-1.7  dative bonds are plain lines; π coordination is drawn to the midpoint of the π bond.

Bond drawing types (bond['t']): 1, 2, 3; 'w' bold wedge; 'h' hashed wedge; 'p' partial double bond;
'd' dashed (partial bond); 'dot' dotted (H bond); 'hp' dashed dative bond; 'bold'; 'arc' chelate backbone sketch.
"""

import math

from ..fonts import adv

D = math.radians


K = dict(lw=0.04 / 0.508, wedge=0.2, sp=0.18, mg=0.056 / 0.508, hs=0.088 / 0.508, fs=0.3528 / 0.508)
LL = 1.5  # 两端都是标签的键


# ---------------- 标签排版 ----------------
# 语法：普通字符；字母或右括号后的数字 → 下标；^x 或 ^{..} → 上标（字母 t/i/n/s 自动斜体）；_x / _{..} → 下标；
#        末尾的 + / − 写成上标电荷。返回字形列表 [(char, x, dy, size, italic)]，以 anchor 字符中心为 0。
def layout(s, fs, anchor=0):
    out = []
    x = 0.0
    i = 0
    k = 0
    prev = ""
    while i < len(s):
        c = s[i]
        if c in "^_" and i + 1 < len(s):
            j = i + 1
            if s[j] == "{":
                e = s.index("}", j)
                seg = s[j + 1 : e]
                i = e + 1
            else:
                seg = s[j]
                i = j + 1
            for ch in seg:
                size = fs * 0.7
                w = adv(ch) * size
                out.append([ch, x, -fs * 0.38 if c == "^" else fs * 0.24, size, c == "^" and ch in "tisn", None])
                x += w
            prev = ""
            continue
        sub = c.isdigit() and prev and (prev.isalpha() or prev in ")]")
        sup = c in "+−" and i == len(s) - 1
        size = fs * (0.7 if (sub or sup) else 1)
        dy = fs * 0.24 if sub else (-fs * 0.42 if sup else 0)
        out.append([c, x, dy, size, False, k])
        x += adv(c) * size
        if not (sub or sup):
            k += 1
        prev = c if not sub else prev
        i += 1
    # anchor：第 anchor 个「正常字符」的中心
    base = [g for g in out if g[5] is not None and g[2] == 0]
    if not base:
        base = out
    a = base[min(anchor, len(base) - 1)]
    cx = a[1] + adv(a[0]) * a[3] / 2
    for g in out:
        g[1] -= cx
    return out


def label_box(atom, B):
    fs = K["fs"] * B
    g = layout(atom["label"], fs, atom.get("anchor", 0))
    x0 = min(q[1] for q in g)
    x1 = max(q[1] + adv(q[0]) * q[3] for q in g)
    m = K["mg"] * B
    y0, y1 = atom["y"] * B - fs * 0.37 - m, atom["y"] * B + fs * 0.37 + m
    # 下标往下多占约 0.24 fs、上标往上多占约 0.42 fs：框要把它们包进去，键才不会压到下标（2001 期 R₂ 的问题）
    if any(q[2] > 0 for q in g):
        y1 += fs * 0.2
    if any(q[2] < 0 for q in g):
        y0 -= fs * 0.15
    if atom.get("hx"):
        if atom.get("hpos", "below") == "below":
            y1 += fs * 1.0
        else:
            y0 -= fs * 1.0
    return (atom["x"] * B + x0 - m, y0, atom["x"] * B + x1 + m, y1)


# ---------------- 分子 ----------------
class M:
    def __init__(s):
        s.at = {}
        s.bd = []
        s.notes = []

    def a(s, k, x, y, label=None, anchor=0, **kw):
        s.at[k] = dict(x=x, y=y, label=label, anchor=anchor, **kw)
        return s

    def go(s, k, frm, ang, label=None, t=1, side=None, L=None, anchor=0, **kw):
        if L is None:
            L = LL if (label and s.at[frm]["label"]) else 1.0
        x0, y0 = s.at[frm]["x"], s.at[frm]["y"]
        s.a(k, x0 + L * math.cos(D(ang)), y0 - L * math.sin(D(ang)), label, anchor, **kw)
        s.b(frm, k, t, side)
        return s

    def b(s, i, j, t=1, side=None, **kw):
        s.bd.append(dict(a=i, b=j, t=t, side=side, **kw))
        return s

    def pos(s, k):
        return (s.at[k]["x"], s.at[k]["y"])

    def ring(s, ids, cx, cy, start, labels=None, side=1.0):
        n = len(ids)
        R = side / (2 * math.sin(math.pi / n))
        labels = labels or {}
        for q, i in enumerate(ids):
            ang = start - 360 / n * q
            if i in s.at:
                continue
            s.a(i, cx + R * math.cos(D(ang)), cy - R * math.sin(D(ang)), labels.get(i))
        return R

    def benzene(s, ids, cx, cy, start, kek=0, side=1.0):
        """苯环：ids 6 个，交替双键（双键画在环内）"""
        s.ring(ids, cx, cy, start, side=side)
        for q in range(6):
            i, j = ids[q], ids[(q + 1) % 6]
            if any({e["a"], e["b"]} == {i, j} for e in s.bd):
                continue
            s.b(i, j, 2 if (q + kek) % 2 == 0 else 1, side="in", ring=ids)
        return s

    def phenyl(s, pre, frm, ang, kek=0):
        """在 frm 上沿 ang 接一个苯基（ipso 到 frm 键长 1，环正六边形）"""
        x0, y0 = s.pos(frm)
        ix, iy = x0 + math.cos(D(ang)), y0 - math.sin(D(ang))
        cx, cy = ix + math.cos(D(ang)), iy - math.sin(D(ang))
        ids = [pre + str(q) for q in range(6)]
        s.a(ids[0], ix, iy)
        s.b(frm, ids[0], 1)
        s.benzene(ids, cx, cy, ang + 180, kek)
        return ids

    def shift(s, dx, dy):
        for a in s.at.values():
            a["x"] += dx
            a["y"] += dy
        return s

    def rot(s, ang, cx=0, cy=0):
        c, n = math.cos(D(ang)), math.sin(D(ang))
        for a in s.at.values():
            x, y = a["x"] - cx, -(a["y"] - cy)
            a["x"], a["y"] = cx + x * c - y * n, cy - (x * n + y * c)
        return s

    def flipx(s):
        for a in s.at.values():
            a["x"] = -a["x"]
        for e in s.bd:
            if isinstance(e.get("side"), (int, float)):
                e["side"] = -e["side"]
        return s


# ---------------- 配位几何模板（IUPAC ST-2） ----------------
# 每个位点：(屏幕角, 键型)。键型 'w' 实楔（指向读者）、'h' 虚楔（背离读者）、1 普通线。
GEOM = {
    # 八面体：轴向上下直线；赤道面四个位点，后方两个虚楔在上（±30°），前方两个实楔在下（夹角 110°）
    "oct": dict(up=(90, 1), dn=(270, 1), bl=(150, "h"), br=(30, "h"), fl=(215, "w"), fr=(325, "w")),
    # 平面四方（透视）：同八面体的赤道面
    "sq": dict(bl=(150, "h"), br=(30, "h"), fl=(215, "w"), fr=(325, "w")),
    # 平面四方（平面画法，四根直线键，与坐标轴成 45°）
    "sqf": dict(ul=(135, 1), ur=(45, 1), dl=(225, 1), dr=(315, 1)),
    # 四面体：两根直线键 + 一实一虚（ST-1.1.1）
    "tet": dict(l=(210, 1), r=(330, 1), w=(105, "w"), h=(75, "h")),
}


def metal(m, k, x, y, sym, geom="oct", rot=0):
    m.a(k, x, y, sym, geom=geom, rot=rot)
    return m


def lig(m, mk, slot, k, label=None, L=None, t=None, anchor=0, **kw):
    g = m.at[mk]
    ang, st = GEOM[g["geom"]][slot]
    ang += g.get("rot", 0)
    if t is not None:
        st = t
    if L is None:
        L = LL if label else 1.0
    m.go(k, mk, ang, label, t=st, L=L, anchor=anchor, **kw)
    return m


# ---------------- 标签锚点 ----------------
def auto_anchor(m):
    """末端标签：键从右边来 → 锚在最后一个正常字符；从左边来或竖直 → 锚在第一个字符（ChemDraw 的做法）。
    已写 anchor_fixed 的原子不动。"""
    nb = {k: [] for k in m.at}
    for e in m.bd:
        nb[e["a"]].append(e["b"])
        nb[e["b"]].append(e["a"])
    for k, a in m.at.items():
        lab = a.get("label")
        if not lab or a.get("anchor_fixed") or len(nb[k]) != 1 or a.get("geom"):
            continue
        o = m.at[nb[k][0]]
        n_norm = len([g for g in layout(lab, 1, 0) if g[5] is not None and g[2] == 0])
        if o["x"] > a["x"] + 0.3:
            a["anchor"] = max(0, n_norm - 1)
        elif o["x"] < a["x"] - 0.3:
            a["anchor"] = 0
    return m


# ---------------- 渲染 ----------------
def f(v):
    return f"{v:.1f}"


def svg_mol(m, B, ox=0, oy=0, ink="#1F1A14", colors=None, bcolors=None):
    colors = colors or {}
    bcolors = bcolors or {}
    auto_anchor(m)
    lw = K["lw"] * B
    out = []
    P = lambda k: (ox + m.at[k]["x"] * B, oy + m.at[k]["y"] * B)
    boxes = {k: (label_box(a, B) if a.get("label") else None) for k, a in m.at.items()}

    def shrink(k, p, q):
        bx = boxes[k]
        if not bx:
            return 0
        bx = (bx[0] + ox, bx[1] + oy, bx[2] + ox, bx[3] + oy)
        dx, dy = q[0] - p[0], q[1] - p[1]
        t = 1e9
        if dx > 1e-9:
            t = min(t, (bx[2] - p[0]) / dx)
        if dx < -1e-9:
            t = min(t, (bx[0] - p[0]) / dx)
        if dy > 1e-9:
            t = min(t, (bx[3] - p[1]) / dy)
        if dy < -1e-9:
            t = min(t, (bx[1] - p[1]) / dy)
        return max(0, min(0.49, t))

    def seg(p, q, i, j):
        t0 = shrink(i, p, q)
        t1 = shrink(j, q, p)
        return (
            (p[0] + (q[0] - p[0]) * t0, p[1] + (q[1] - p[1]) * t0),
            (q[0] + (p[0] - q[0]) * t1, q[1] + (p[1] - q[1]) * t1),
        )

    line = lambda p, q, c, w=lw, da="": (
        f'<line x1="{f(p[0])}" y1="{f(p[1])}" x2="{f(q[0])}" y2="{f(q[1])}" stroke="{c}" stroke-width="{f(w)}" stroke-linecap="round"{(" stroke-dasharray=%s" % chr(34) + da + chr(34)) if da else ""}/>'
    )
    for e in sorted(m.bd, key=lambda e: e.get("z", 0)):
        i, j, t = e["a"], e["b"], e["t"]
        if e.get("halo"):  # 前面的键压在后面的键上：先铺一条白色宽线，形成 IUPAC ST-3 的「断开」效果
            pa, pb = P(i), P(j)
            a0, b0 = seg(pa, pb, i, j)
            hw = (K["sp"] * B if t in (2, 3, "p") else 0) + K["lw"] * B * 5
            ux_, uy_ = b0[0] - a0[0], b0[1] - a0[1]
            LL_ = math.hypot(ux_, uy_) or 1
            cut = K["lw"] * B * 3
            a1 = (a0[0] + ux_ / LL_ * cut, a0[1] + uy_ / LL_ * cut)
            b1 = (b0[0] - ux_ / LL_ * cut, b0[1] - uy_ / LL_ * cut)
            out.append(
                f'<line x1="{f(a1[0])}" y1="{f(a1[1])}" x2="{f(b1[0])}" y2="{f(b1[1])}" stroke="#FFFFFF" stroke-width="{f(hw)}"/>'
            )
        c = bcolors.get((i, j)) or bcolors.get((j, i)) or e.get("col") or ink
        p, q = P(i), P(j)
        dx, dy = q[0] - p[0], q[1] - p[1]
        L = math.hypot(dx, dy)
        ux, uy = dx / L, dy / L
        nx, ny = -uy, ux
        a, b = seg(p, q, i, j)
        if t == "w":
            w = K["wedge"] * B / 2
            w0 = lw * 0.45
            out.append(
                f'<path d="M{f(a[0] + nx * w0)} {f(a[1] + ny * w0)} L{f(b[0] + nx * w)} {f(b[1] + ny * w)} L{f(b[0] - nx * w)} {f(b[1] - ny * w)} L{f(a[0] - nx * w0)} {f(a[1] - ny * w0)}Z" fill="{c}" stroke="{c}" stroke-width="{f(lw * 0.5)}" stroke-linejoin="round"/>'
            )
            continue
        if t == "h":
            Lq = math.hypot(b[0] - a[0], b[1] - a[1])
            n = max(3, round(Lq / (K["hs"] * B)))
            for s in range(1, n + 1):
                u = s / n
                cx, cy = a[0] + (b[0] - a[0]) * u, a[1] + (b[1] - a[1]) * u
                w = max(lw * 0.6, K["wedge"] * B * u) / 2
                out.append(line((cx + nx * w, cy + ny * w), (cx - nx * w, cy - ny * w), c, lw * 0.85))
            continue
        if t == "hp":  # 不带楔形的虚线配位（等宽短线），用于背后的配位
            Lq = math.hypot(b[0] - a[0], b[1] - a[1])
            n = max(3, round(Lq / (K["hs"] * B)))
            for s in range(0, n + 1):
                u = s / n
                cx, cy = a[0] + (b[0] - a[0]) * u, a[1] + (b[1] - a[1]) * u
                w = K["wedge"] * B * 0.32
                out.append(line((cx + nx * w, cy + ny * w), (cx - nx * w, cy - ny * w), c, lw * 0.85))
            continue
        if t == "arc":  # 螯合骨架的弧线（如 BINAP 的联萘骨架示意），bulge 为向外鼓的键长倍数，正 = 法向左侧
            bg = e.get("bulge", 0.6) * B
            cx_, cy_ = (p[0] + q[0]) / 2 + nx * bg, (p[1] + q[1]) / 2 + ny * bg

            def sh(P0, k):
                bx = boxes[k]
                if not bx:
                    return P0
                vx, vy = cx_ - P0[0], cy_ - P0[1]
                l = math.hypot(vx, vy)
                r = K["fs"] * B * 0.5
                return (P0[0] + vx / l * r, P0[1] + vy / l * r)

            A0 = sh(p, i)
            B0 = sh(q, j)
            out.append(
                f'<path d="M{f(A0[0])} {f(A0[1])} Q{f(cx_)} {f(cy_)} {f(B0[0])} {f(B0[1])}" fill="none" stroke="{c}" stroke-width="{f(lw)}"/>'
            )
            continue
        if t == "bold":
            out.append(line(a, b, c, K["wedge"] * B * 0.9))
            continue
        if t == "d":  # 部分键 / 氢键：虚线
            out.append(line(a, b, c, lw, f"{f(0.10 * B)} {f(0.08 * B)}"))
            continue
        if t == "dot":  # 氢键：圆点
            Lq = math.hypot(b[0] - a[0], b[1] - a[1])
            n = max(3, int(Lq / (0.13 * B)))
            for s_ in range(n + 1):
                u = s_ / n
                out.append(
                    f'<circle cx="{f(a[0] + (b[0] - a[0]) * u)}" cy="{f(a[1] + (b[1] - a[1]) * u)}" r="{f(lw * 0.75)}" fill="{c}"/>'
                )
            continue
        if t == 1:
            out.append(line(a, b, c))
            continue
        sp = K["sp"] * B
        side = e.get("side")
        # 连着末端原子（=O、=NH、=CH₂ 等只有这一根键）的双键一律居中画（ChemDraw 同样处理），两线等长平齐；
        # 偏向一侧只留给环内双键和链中间的 C=C（用户 2026-10-16：双键不要一高一低）
        if side not in (0, None, "in") and t in (2, "p") and not e.get("ring"):
            deg = lambda k: sum(1 for e2 in m.bd if k in (e2["a"], e2["b"]))
            if deg(i) == 1 or deg(j) == 1:
                side = None
        if side == "in":  # 环内双键：朝环心
            ids = e["ring"]
            gx = sum(m.at[k]["x"] for k in ids) / len(ids) * B + ox
            gy = sum(m.at[k]["y"] for k in ids) / len(ids) * B + oy
            mx, my = (p[0] + q[0]) / 2, (p[1] + q[1]) / 2
            side = 1 if (gx - mx) * nx + (gy - my) * ny > 0 else -1

        # 多重键的几条平行线共用同一对截断参数（用户 2026-10-16：「双键总是一高一低……在任何倾斜角度下，两个键都要同样平齐，不要错开」）。
        # 原来每条线各自被标签矩形截断，斜着的时候靠矩形角的那条截得多、另一条截得少，两端就错开了。
        # 现在：每端取所有平行线里最大的截断量（保证都不压标签），所有线都在同一个 [t0, t1] 上画，端点连线垂直于键轴。
        def par(offs):
            P2 = [((p[0] + nx * sp * o, p[1] + ny * sp * o), (q[0] + nx * sp * o, q[1] + ny * sp * o)) for o in offs]
            t0 = max(shrink(i, A_, B_) for A_, B_ in P2)
            t1 = max(shrink(j, B_, A_) for A_, B_ in P2)
            return [((A_[0] + dx * t0, A_[1] + dy * t0), (B_[0] - dx * t1, B_[1] - dy * t1)) for A_, B_ in P2]

        if t == 3:
            for a2, b2 in par((-1, 0, 1)):
                out.append(line(a2, b2, c))
            continue
        if side in (0, None):
            for (a2, b2), sd in zip(par((-0.5, 0.5)), (-0.5, 0.5)):
                out.append(line(a2, b2, c, lw, f"{f(0.09 * B)} {f(0.07 * B)}" if (t == "p" and sd > 0) else ""))
            continue
        (a, b), (a2, b2) = par((0, side))
        out.append(line(a, b, c))
        # 偏向一侧的内线（环内双键）：在与外线平齐的基础上两端等量缩进，左右对称，不出现一端齐、一端缩的错位
        sh = sp * 0.85
        a2 = (a2[0] + ux * sh, a2[1] + uy * sh)
        b2 = (b2[0] - ux * sh, b2[1] - uy * sh)
        out.append(line(a2, b2, c, lw, f"{f(0.09 * B)} {f(0.07 * B)}" if t == "p" else ""))
    fs = K["fs"] * B
    for k, a in m.at.items():
        if a.get("label") and a.get("hx"):
            col = colors.get(k) or a.get("col") or ink
            g = layout(a["hx"], fs, 0)
            x = ox + a["x"] * B - adv(a["hx"][0]) * fs / 2
            y = oy + a["y"] * B + fs * 0.36 + (fs * 1.0 if a.get("hpos", "below") == "below" else -fs * 1.0)
            sp = "".join(
                f'<tspan x="{f(x + q[1] + adv(a["hx"][0]) * fs / 2)}" y="{f(y + q[2])}" font-size="{f(q[3])}">{q[0]}</tspan>'
                for q in g
            )
            out.append(f'<text fill="{col}" style="font-family:\'Liberation Sans\',Arial,sans-serif">{sp}</text>')
    for k, a in m.at.items():
        if a.get("label") and a.get("chg"):  # 电荷放在标签左上（避开右边的键）
            out.append(
                f'<text x="{f(ox + a["x"] * B - fs * 0.62)}" y="{f(oy + a["y"] * B - fs * 0.2)}" font-size="{f(fs * 0.7)}" fill="{ink}" text-anchor="middle" style="font-family:\'Liberation Sans\',Arial,sans-serif">{a["chg"]}</text>'
            )
    for k, a in m.at.items():
        if not a.get("label"):
            continue
        col = colors.get(k) or a.get("col") or ink
        g = layout(a["label"], fs, a.get("anchor", 0))
        x, y = ox + a["x"] * B, oy + a["y"] * B + fs * 0.36
        loc = locant_mask(g)
        sp = "".join(
            f'<tspan x="{f(x + q[1])}" y="{f(y + q[2])}" font-size="{f(q[3])}"{" font-style=%sitalic%s" % (chr(34), chr(34)) if q[4] else ""}{(" fill=%s%s%s" % (chr(34), LOC_COL, chr(34))) if loc[ii] else ""}>{"&lt;" if q[0] == "<" else q[0]}</tspan>'
            for ii, q in enumerate(g)
        )
        out.append(f'<text fill="{col}" style="font-family:\'Liberation Sans\',Arial,sans-serif">{sp}</text>')
    return "\n".join(out)


# v2（2026-10-14 用户：「对碳原子的顺序标识数字可以用另一种颜色」）：标签里纯数字的上标（^{6}CH2、HO^{3}CH）是原子编号，
# 设了 LOC_COL 就单独上色；带正负号的上标（^{2−}）是电荷，不动。LOC_COL = None 时与 v1 完全一样。
LOC_COL = None


def locant_mask(g):
    m = [False] * len(g)
    if not LOC_COL:
        return m
    i = 0
    while i < len(g):
        if g[i][5] is None and g[i][2] < 0:
            j = i
            while j < len(g) and g[j][5] is None and g[j][2] < 0:
                j += 1
            run = "".join(q[0] for q in g[i:j])
            nxt = g[j][0] if j < len(g) else ""
            prv = g[i - 1][0] if i > 0 else ""
            # 原子编号的上标后面紧跟原子符号（^{6}CH2、HO^{3}CH）；残基号 Asp^{25}、R^{1} 后面不是大写字母 → 不上色（r4/hiv 组）
            if (
                run
                and all(c.isdigit() or c in "′',α-βγ" for c in run)
                and any(c.isdigit() for c in run)
                and nxt[:1].isupper()
                and prv in ("", "H", "O", "N", "S", ")")
            ):
                for t in range(i, j):
                    m[t] = True
            i = j
        else:
            i += 1
    return m


# ---------------- 几何检查 ----------------
def check(m, name="", B=40, allow=()):
    """返回问题列表：键长、键角（只查普通键）、楔形窄端、标签与键 / 标签互相压住、键交叉"""
    issues = []
    auto_anchor(m)
    P = lambda k: (m.at[k]["x"], m.at[k]["y"])
    nb = {k: [] for k in m.at}
    for e in m.bd:
        if e["t"] == "arc":
            continue
        nb[e["a"]].append(e)
        nb[e["b"]].append(e)
    # 键长
    for e in m.bd:
        if e["t"] in ("d", "dot", "hp", "arc") or e.get("free") or m.at[e["a"]].get("geom") or m.at[e["b"]].get("geom"):
            continue
        L = math.dist(P(e["a"]), P(e["b"]))
        la, lb = m.at[e["a"]].get("label"), m.at[e["b"]].get("label")
        want = (1.0, 1.2, 1.3, LL) if (la and lb) else (1.0,)
        if "L0" in e:
            want = (e["L0"],)
        if min(abs(L - w) for w in want) > 0.04 and not e.get("len_ok"):
            issues.append(f"{name}: 键 {e['a']}–{e['b']} 长 {L:.2f}")
    # 键角：对每个非金属原子，所有画出来的键（含楔形）两两相邻夹角
    metals = {k for k, a in m.at.items() if a.get("geom")}
    for k, es in nb.items():
        if k in metals or k in allow:
            continue
        if len(es) < 2:
            continue
        angs = []
        for e in es:
            o = e["b"] if e["a"] == k else e["a"]
            if e["t"] in ("d", "dot") or (o in metals and e["t"] == "hp"):
                continue
            x0, y0 = P(k)
            x1, y1 = P(o)
            angs.append((math.degrees(math.atan2(-(y1 - y0), x1 - x0)) % 360, e["t"], o))
        angs.sort()
        n = len(angs)
        if n < 2:
            continue
        gaps = [((angs[(q + 1) % n][0] - angs[q][0]) % 360, angs[q][2], angs[(q + 1) % n][2]) for q in range(n)]
        if n == 2:
            g = min(gaps[0][0], 360 - gaps[0][0])
            ts = [e["t"] for e in es]
            sp = 3 in ts or ts.count(2) == 2 or m.at[k].get("linear")  # 只有 sp 原子（连三键 / 累积双键）才画 180°
            ok = (g >= 174) if sp else (112 <= g <= 128)
            if not (ok or m.at[k].get("ring_ok")):
                issues.append(f"{name}: {k} 两键夹角 {g:.0f}°（{'sp 应 180°' if sp else '应 120°'}）")
        elif n == 3:
            if not m.at[k].get("ring_ok"):
                for g, a1, a2 in gaps:
                    if not (108 <= g <= 132):
                        issues.append(f"{name}: {k} 三键夹角 {a1}–{k}–{a2} = {g:.0f}°")
        else:
            for g, a1, a2 in gaps:
                if g < 45:
                    issues.append(f"{name}: {k} 键夹角过小 {a1}–{k}–{a2} = {g:.0f}°")
    # 看得见的键长：标签留白剪掉之后，画出来的那段 ≥ 0.4 键长（否则是“短桩”，如 O–O、Os–O）
    for e in m.bd:
        if e["t"] in ("arc", "dot") or m.at[e["a"]].get("pi") or m.at[e["b"]].get("pi"):
            continue
        pa, pb = P(e["a"]), P(e["b"])
        L = math.dist(pa, pb)
        if L < 1e-6:
            continue

        def _cut(k, p, q):
            a = m.at[k]
            if not a.get("label"):
                return 0
            bx = label_box(a, 1)
            dx, dy = q[0] - p[0], q[1] - p[1]
            t = 1e9
            if dx > 1e-9:
                t = min(t, (bx[2] - p[0]) / dx)
            if dx < -1e-9:
                t = min(t, (bx[0] - p[0]) / dx)
            if dy > 1e-9:
                t = min(t, (bx[3] - p[1]) / dy)
            if dy < -1e-9:
                t = min(t, (bx[1] - p[1]) / dy)
            return max(0, min(0.49, t))

        vis = L * (1 - _cut(e["a"], pa, pb) - _cut(e["b"], pb, pa))
        metal_b = m.at[e["a"]].get("geom") or m.at[e["b"]].get("geom")
        both = m.at[e["a"]].get("label") and m.at[e["b"]].get("label")
        lim = (
            0.33 if metal_b else (0.425 if both else 0.4)
        )  # 第四轮：两端都是标签的键可见段 ≥ 0.43（plp 评审：0.42 的 O5′ 短桩挂小钩）
        if vis < lim and not e.get("stub_ok"):
            issues.append(f"{name}: 键 {e['a']}–{e['b']} 画出来只有 {vis:.2f}（短桩，应 ≥ {lim}）")
    # 楔形窄端：wedge/hash 的 a 端应是立体中心（有 ≥3 个邻居）或金属
    for e in m.bd:
        if e["t"] in ("w", "h"):
            if len(nb[e["a"]]) < 3 and e["a"] not in metals:
                issues.append(f"{name}: 楔形 {e['a']}→{e['b']} 窄端不在立体中心")
    # 标签框与无关的键相交、标签互相重叠
    boxes = {k: label_box(a, 1) for k, a in m.at.items() if a.get("label")}

    def seg_box(p, q, bx):
        # 线段与矩形是否相交（采样）
        for s in range(0, 21):
            u = s / 20
            x = p[0] + (q[0] - p[0]) * u
            y = p[1] + (q[1] - p[1]) * u
            if bx[0] < x < bx[2] and bx[1] < y < bx[3]:
                return True
        return False

    for k, bx in boxes.items():
        shr = (bx[0] + 0.03, bx[1] + 0.03, bx[2] - 0.03, bx[3] - 0.03)
        for e in m.bd:
            if k in (e["a"], e["b"]) or e["t"] == "arc":
                continue
            if seg_box(P(e["a"]), P(e["b"]), shr):
                issues.append(f"{name}: 标签 {k}({m.at[k]['label']}) 压住键 {e['a']}–{e['b']}")
    ks = list(boxes)
    for q in range(len(ks)):
        for r in range(q + 1, len(ks)):
            a, b = boxes[ks[q]], boxes[ks[r]]
            if a[0] < b[2] - 0.02 and b[0] < a[2] - 0.02 and a[1] < b[3] - 0.02 and b[1] < a[3] - 0.02:
                issues.append(f"{name}: 标签 {ks[q]} 与 {ks[r]} 重叠")

    # 键交叉
    def inter(p1, p2, p3, p4):
        d = lambda a, b, c: (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
        return d(p1, p2, p3) * d(p1, p2, p4) < 0 and d(p3, p4, p1) * d(p3, p4, p2) < 0

    for q in range(len(m.bd)):
        for r in range(q + 1, len(m.bd)):
            e1, e2 = m.bd[q], m.bd[r]
            if {e1["a"], e1["b"]} & {e2["a"], e2["b"]}:
                continue
            if e1["t"] == "arc" or e2["t"] == "arc":
                continue
            # π 配位线落在 C=C / C=N 中点上，本来就“碰”那根键
            pis = [m.at[x].get("pi") for x in (e1["a"], e1["b"], e2["a"], e2["b"]) if m.at[x].get("pi")]
            if any(set(pp) in ({e1["a"], e1["b"]}, {e2["a"], e2["b"]}) for pp in pis):
                continue
            if inter(P(e1["a"]), P(e1["b"]), P(e2["a"]), P(e2["b"])):
                if not (e1.get("cross_ok") or e2.get("cross_ok") or e1.get("halo") or e2.get("halo")):
                    issues.append(f"{name}: 键 {e1['a']}–{e1['b']} 与 {e2['a']}–{e2['b']} 交叉")
    return issues


def bbox(m, B=1):
    auto_anchor(m)
    xs, ys = [], []
    for k, a in m.at.items():
        if a.get("label"):
            bx = label_box(a, 1)
            xs += [bx[0], bx[2]]
            ys += [bx[1], bx[3]]
        else:
            xs.append(a["x"])
            ys.append(a["y"])
    return min(xs) * B, min(ys) * B, max(xs) * B, max(ys) * B


# ---------------- 2D 松弛（把粗排的坐标整理成标准键长 / 键角，同时推开重叠的标签） ----------------
def relax(m, fixed=(), anchor_w=0.02, iters=3000, coord=("hp", "d", "dot"), ang_w=1.0, verbose=False):
    """fixed：不动的原子。键长目标：bond['L0'] 或 1.0（金属键与 coord 类型键不约束长度，除非给了 L0）。
    键角目标：每个非金属原子的「有机键」两两相邻夹角 → 2 键 120°（atom['lin'] 时 180°），3 键 120°；atom['ang_off'] 关掉。
    另加：标签互不重叠、标签不压无关的键、非键原子间距 ≥ 0.75。"""
    import numpy as np
    from scipy.optimize import minimize

    ks = list(m.at)
    idx = {k: i for i, k in enumerate(ks)}
    metals = {k for k, a in m.at.items() if a.get("geom")}
    pis = {k: a["pi"] for k, a in m.at.items() if a.get("pi")}
    free = [k for k in ks if k not in fixed and k not in metals and k not in pis]
    fi = {k: i for i, k in enumerate(free)}
    X0 = np.array([[m.at[k]["x"], m.at[k]["y"]] for k in ks])
    # 键长项
    LB = []
    CR = []
    for e in m.bd:
        a, b = e["a"], e["b"]
        if "L0" in e:
            LB.append((idx[a], idx[b], e["L0"], e.get("Lw", 4.0)))
        elif (a in metals or b in metals) and not e.get("fixedlen"):
            if a in fixed and b in fixed:
                continue
            CR.append((idx[a], idx[b], e.get("Lmin", 1.5), e.get("Lmax", 2.1)))
        elif e["t"] in coord:
            continue
        else:
            LB.append((idx[a], idx[b], e.get("L0", 1.0), 40.0))
    # 键角项
    nb = {k: [] for k in ks}
    for e in m.bd:
        if e["t"] in ("d", "dot", "arc"):
            continue
        if (e["a"] in metals or e["b"] in metals) and not e.get("ang"):
            continue
        nb[e["a"]].append(e["b"])
        nb[e["b"]].append(e["a"])
    AT = []
    for k, ns in nb.items():
        a = m.at[k]
        if k in metals or a.get("ang_off") or len(ns) < 2:
            continue
        if len(ns) == 2:
            AT.append((idx[k], [idx[n] for n in ns], 180 if a.get("lin") else a.get("ang2", 120)))
        elif len(ns) == 3:
            AT.append((idx[k], [idx[n] for n in ns], 120))
    lab = [k for k in ks if m.at[k].get("label")]
    boxes0 = {k: label_box(dict(m.at[k], x=0, y=0), 1) for k in lab}
    bonds = [(idx[e["a"]], idx[e["b"]]) for e in m.bd if e["t"] != "arc"]

    def unpack(v):
        X = X0.copy()
        for k, i in fi.items():
            X[idx[k]] = v[2 * i : 2 * i + 2]
        for k, (p, q) in pis.items():
            X[idx[k]] = (X[idx[p]] + X[idx[q]]) / 2
        return X

    # 预先整理成数组，能量函数全部向量化
    LBi = np.array([[a, b] for a, b, _, _ in LB], int).reshape(-1, 2)
    LBL = np.array([l for _, _, l, _ in LB])
    LBw = np.array([w for *_, w in LB])
    A2 = [(c, ns, t) for c, ns, t in AT if len(ns) == 2]
    A3 = [(c, ns) for c, ns, t in AT if len(ns) == 3]
    A2c = np.array([c for c, _, _ in A2], int)
    A2n = np.array([ns for _, ns, _ in A2], int).reshape(-1, 2)
    A2t = np.radians([t for *_, t in A2])
    A3c = np.array([c for c, _ in A3], int)
    A3n = np.array([ns for _, ns in A3], int).reshape(-1, 3)
    Li = np.array([idx[k] for k in lab], int)
    Bx = np.array([boxes0[k] for k in lab]).reshape(-1, 4)
    BD = np.array(bonds, int).reshape(-1, 2)
    us = np.array([0.15, 0.3, 0.5, 0.7, 0.85])
    inc = np.array([[i0 in (i, j) for i, j in bonds] for i0 in Li], bool).reshape(len(Li), len(bonds))
    n = len(ks)
    bonded = np.zeros((n, n), bool)
    for i, j in bonds:
        bonded[i, j] = bonded[j, i] = True
    pis_i = [idx[k] for k in pis]
    nb_mask = ~bonded & np.triu(np.ones((n, n), bool), 1)
    for p in pis_i:
        nb_mask[p, :] = False
        nb_mask[:, p] = False
    pair_q, pair_r = np.where(nb_mask)
    lq, lr = np.triu_indices(len(Li), 1)
    fidx = np.array([idx[k] for k in free], int)
    CRi = np.array([[a, b] for a, b, _, _ in CR], int).reshape(-1, 2)
    CRlo = np.array([l for _, _, l, _ in CR])
    CRhi = np.array([h for *_, h in CR])
    # 键—键接近：非相邻、非 halo 的键对，取样点距离 < 0.3 罚
    BB = []
    bd2 = [e for e in m.bd if e["t"] != "arc"]
    for q in range(len(bd2)):
        for r in range(q + 1, len(bd2)):
            e1, e2 = bd2[q], bd2[r]
            if {e1["a"], e1["b"]} & {e2["a"], e2["b"]}:
                continue
            if e1.get("halo") or e2.get("halo") or e1.get("cross_ok") or e2.get("cross_ok"):
                continue
            if m.at[e1["a"]].get("pi") or m.at[e1["b"]].get("pi") or m.at[e2["a"]].get("pi") or m.at[e2["b"]].get("pi"):
                pass
            BB.append((q, r))
    BBq = np.array([q for q, _ in BB], int)
    BBr = np.array([r for _, r in BB], int)
    us2 = np.linspace(0.1, 0.9, 5)

    def E(v):
        X = unpack(v)
        e = 0.0
        if len(LBi):
            d = np.linalg.norm(X[LBi[:, 0]] - X[LBi[:, 1]], axis=1)
            e += np.sum(LBw * (d - LBL) ** 2)
        if len(A2c):
            u1 = X[A2n[:, 0]] - X[A2c]
            u2 = X[A2n[:, 1]] - X[A2c]
            cos = np.sum(u1 * u2, 1) / (np.linalg.norm(u1, axis=1) * np.linalg.norm(u2, axis=1) + 1e-9)
            g = np.arccos(np.clip(cos, -1, 1))
            e += 2 * ang_w * np.sum(((g - A2t) / np.radians(30)) ** 2)
        if len(A3c):
            ang = np.arctan2(-(X[A3n][:, :, 1] - X[A3c][:, None, 1]), X[A3n][:, :, 0] - X[A3c][:, None, 0])
            ang = np.sort(ang % (2 * np.pi), 1)
            gaps = np.diff(np.concatenate([ang, ang[:, :1] + 2 * np.pi], 1), axis=1)
            e += 2 * ang_w * np.sum(((gaps - np.radians(120)) / np.radians(30)) ** 2)
        if len(Li) > 1:
            P = X[Li]
            lo = P + Bx[:, :2]
            hi = P + Bx[:, 2:]
            ox = np.minimum(hi[lq, 0], hi[lr, 0]) - np.maximum(lo[lq, 0], lo[lr, 0])
            oy = np.minimum(hi[lq, 1], hi[lr, 1]) - np.maximum(lo[lq, 1], lo[lr, 1])
            mk = (ox > -0.18) & (oy > -0.18)
            e += 60 * np.sum((ox[mk] + 0.18) * (oy[mk] + 0.18))
        if len(Li) and len(BD):
            P = X[Li]
            c = P + (Bx[:, :2] + Bx[:, 2:]) / 2
            r = (Bx[:, 2:] - Bx[:, :2]) / 2 + 0.05
            Q = X[BD[:, 0]][:, None, :] + (X[BD[:, 1]] - X[BD[:, 0]])[:, None, :] * us[None, :, None]  # bonds × 5 × 2
            dx = np.abs(Q[None, :, :, 0] - c[:, None, None, 0]) / r[:, None, None, 0]
            dy = np.abs(Q[None, :, :, 1] - c[:, None, None, 1]) / r[:, None, None, 1]
            mm = np.maximum(dx, dy)
            pen = np.where(mm < 1, (1 - mm) ** 2, 0) * (~inc)[:, :, None]
            e += 25 * np.sum(pen)
        if len(pair_q):
            d = np.linalg.norm(X[pair_q] - X[pair_r], axis=1)
            e += 5 * np.sum(np.clip(0.75 - d, 0, None) ** 2)
        if len(CRi):
            d = np.linalg.norm(X[CRi[:, 0]] - X[CRi[:, 1]], axis=1)
            e += 10 * np.sum(np.clip(CRlo - d, 0, None) ** 2) + 4 * np.sum(np.clip(d - CRhi, 0, None) ** 2)
        if len(BBq):
            S0 = X[BD[:, 0]]
            S1 = X[BD[:, 1]]
            Pq = S0[BBq][:, None, :] + (S1[BBq] - S0[BBq])[:, None, :] * us2[None, :, None]
            Pr = S0[BBr][:, None, :] + (S1[BBr] - S0[BBr])[:, None, :] * us2[None, :, None]
            dd = np.linalg.norm(Pq[:, :, None, :] - Pr[:, None, :, :], axis=3)
            e += 30 * np.sum(np.clip(0.32 - dd, 0, None) ** 2)
        e += anchor_w * np.sum((X[fidx] - X0[fidx]) ** 2)
        return e

    v0 = np.concatenate([X0[idx[k]] for k in free]) if free else np.zeros(0)
    if not free:
        return m
    r = minimize(E, v0, method="L-BFGS-B", options=dict(maxiter=iters))
    X = unpack(r.x)
    for k in free + list(pis):
        m.at[k]["x"], m.at[k]["y"] = float(X[idx[k]][0]), float(X[idx[k]][1])
    m.relaxE = float(r.fun)
    if verbose:
        print("relax E", E(v0), "→", r.fun)
    return m
