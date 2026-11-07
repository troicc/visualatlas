"""Read the CIP label of a stereocentre back from a drawing (2D coordinates + wedges).

The SMILES uses atom-map numbers ([C:1]) to say which drawn atom each SMILES atom is; atoms that are not
drawn (condensed into a label) are placed near the atom they are attached to."""

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdCIPLabeler
from rdkit.Geometry import Point3D


def cip_from_drawing(m, smiles, keymap, center_key, replace=None):
    """keymap：{映射号: 图中键名}；replace：图中键名 → 在产物里代替它的原子（如 Rh → H，用映射号指明位置）"""
    _p = Chem.SmilesParserParams()
    _p.removeHs = False
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles, _p))
    n = mol.GetNumAtoms()
    pos = {}
    inv = {}
    for a in mol.GetAtoms():
        mp = a.GetAtomMapNum()
        if mp and mp in keymap:
            inv[a.GetIdx()] = keymap[mp]
    # 中心与邻居的 z
    ck = center_key
    zmap = {ck: 0.0}
    for e in m.bd:
        if ck in (e["a"], e["b"]):
            o = e["b"] if e["a"] == ck else e["a"]
            t = e["t"]
            if e["a"] == ck:
                z = 0.8 if t == "w" else (-0.8 if t == "h" else 0)
            else:
                z = -0.8 if t == "w" else (0.8 if t == "h" else 0)  # 窄端在对方：方向相反
            zmap[o] = z
    for i, k in inv.items():
        a = m.at[k]
        pos[i] = np.array([a["x"], -a["y"], zmap.get(k, 0.0)])
    # 守卫：中心的重原子邻居必须都映射到图中，否则结果无意义
    for a in mol.GetAtoms():
        if a.GetAtomMapNum() and keymap.get(a.GetAtomMapNum()) == ck:
            for x in a.GetNeighbors():
                if x.GetAtomicNum() > 1 and not (x.GetAtomMapNum() and x.GetAtomMapNum() in keymap):
                    raise ValueError(f"立体中心 {ck} 的邻居 {x.GetSymbol()}{x.GetIdx()} 没有映射")
    # 中心上没画出来的 H：放在已画邻居向量和的反方向
    for a in mol.GetAtoms():
        if a.GetAtomMapNum() and keymap.get(a.GetAtomMapNum()) == ck:
            c = a.GetIdx()
            vs = [pos[x.GetIdx()] - pos[c] for x in a.GetNeighbors() if x.GetIdx() in pos]
            for x in a.GetNeighbors():
                if x.GetIdx() not in pos and x.GetAtomicNum() == 1:
                    v = -sum(vs)
                    pos[x.GetIdx()] = pos[c] + v / (np.linalg.norm(v) + 1e-9)
    # 其余原子：放在已定位邻居附近（随机小偏移）
    rng = np.random.default_rng(0)
    for _ in range(n):
        for a in mol.GetAtoms():
            i = a.GetIdx()
            if i in pos:
                continue
            nb = [x.GetIdx() for x in a.GetNeighbors() if x.GetIdx() in pos]
            if nb:
                pos[i] = pos[nb[0]] + rng.normal(0, 0.3, 3)
    conf = Chem.Conformer(n)
    for i in range(n):
        conf.SetAtomPosition(i, Point3D(*map(float, pos.get(i, np.zeros(3)))))
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    Chem.AssignStereochemistryFrom3D(mol, replaceExistingTags=True)
    rdCIPLabeler.AssignCIPLabels(mol)
    for a in mol.GetAtoms():
        if a.GetAtomMapNum() and keymap.get(a.GetAtomMapNum()) == ck:
            return a.GetProp("_CIPCode") if a.HasProp("_CIPCode") else None
