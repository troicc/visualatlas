"""Baseline 1 (PLAN.md, experiment 1): render mechanisms with ChRIMP's own visualizer.

Runs in ChRIMP's environment, not ours (see baselines/README.md):

    python chrimp_render.py <ChRIMP checkout> <out dir> [rxn_idx ...]

Draws every listed reaction of ChRIMP's preprocessed mech-USPTO-31K as one multi-step figure
(MechanismVisualizer) and each step on its own (MechSmiles.show_reac).
"""
import contextlib
import io
import os
import sys

import pandas as pd
from chrimp.visualization.mechanism_visualizer import MechanismVisualizer
from chrimp.world.mechsmiles import MechSmiles

root, out = sys.argv[1], sys.argv[2]
rxns = [int(x) for x in sys.argv[3:]] or [0, 3, 5]
os.makedirs(out, exist_ok=True)
df = pd.read_csv(os.path.join(root, "data", "preprocessed", "mech_uspto_31k.csv"))

for r in rxns:
    # step_idx_retro counts backwards from the product
    steps = df[df.rxn_idx == r].sort_values("step_idx_retro", ascending=False).mech_smi.tolist()
    with contextlib.redirect_stdout(io.StringIO()):
        MechanismVisualizer(steps).show(save_path=os.path.join(out, f"rxn{r}.svg"), return_svg=True)
        for k, s in enumerate(steps):
            svg = MechSmiles(s).show_reac(return_svg=True)
            with open(os.path.join(out, f"rxn{r}_step{k}.svg"), "w") as fh:
                fh.write(svg if isinstance(svg, str) else svg.decode())
    print(f"rxn {r}: {len(steps)} steps")
