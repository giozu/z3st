#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Damage field at three steps, and an animation of the crack growth.

Reads `output/fields_*.vtu` and `output/force_displacement.csv` (run
`plot_force_displacement.py` first) and writes `output/damage_steps.png` and
`output/damage_evolution.gif`. The case-specific values are the constants
below.
"""

import os
from glob import glob

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import pyvista as pv
from matplotlib.animation import FuncAnimation, PillowWriter

N_STEPS = 701
STEPS = [550, 625, 700]          # panels of damage_steps.png
FRAMES = range(450, 701, 3)      # before step 450 the damage stays on the notch
U_LABEL = "u_y"                  # loading direction

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
VTU_FILES = sorted(glob(os.path.join(OUTPUT_DIR, "fields_*.vtu")))
if len(VTU_FILES) < N_STEPS:
    raise FileNotFoundError(f"Need the {N_STEPS} VTU files of a full run in {OUTPUT_DIR}")
u = np.loadtxt(os.path.join(OUTPUT_DIR, "force_displacement.csv"),
               delimiter=",", skiprows=1)[:, 1] * 1e3  # um

m0 = pv.read(VTU_FILES[0])
# Triangles from the mesh cells: a Delaunay triangulation would bridge the slit.
tri = mtri.Triangulation(m0.points[:, 0] * 1e3, m0.points[:, 1] * 1e3,
                         np.asarray(m0.cell_connectivity).reshape(-1, 3))
levels = np.linspace(0.0, 1.0, 21)


def draw(ax, i):
    d = np.asarray(pv.read(VTU_FILES[i]).point_data["Damage"]).ravel()
    cf = ax.tricontourf(tri, d, levels=levels, cmap="hot_r")
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title(f"step {i}, ${U_LABEL}$ = {u[i]:.2f} µm")
    return cf


fig, axs = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
for ax, i in zip(axs, STEPS):
    cf = draw(ax, i)
fig.colorbar(cf, ax=axs, label="Damage $d$", shrink=0.8)
fig.savefig(os.path.join(OUTPUT_DIR, "damage_steps.png"), dpi=200)
plt.close(fig)
print("[INFO] damage_steps.png saved")

fig, ax = plt.subplots(figsize=(6, 5.5), constrained_layout=True)
fig.colorbar(draw(ax, 0), ax=ax, label="Damage $d$")


def update(i):
    for c in list(ax.collections):
        c.remove()
    draw(ax, i)


frames = list(FRAMES) + [N_STEPS - 1]
FuncAnimation(fig, update, frames=frames).save(
    os.path.join(OUTPUT_DIR, "damage_evolution.gif"), writer=PillowWriter(fps=12), dpi=80)
print("[INFO] damage_evolution.gif saved")
