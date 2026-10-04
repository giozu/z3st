#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""Figures of the plate-with-hole case, written to output/.

  mesh.png           the quarter plate, and a close-up of the refined mesh at the hole
  stress_field.png   sigma_xx and von Mises around the hole, within 3 hole radii
  kirsch.png         sigma_xx across the section x = 0 and the hoop stress on the
                     hole edge, Z3ST against Kirsch

Run after the case: python3 plots.py (Allrun does it).
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import meshio  # noqa: E402
import numpy as np  # noqa: E402
import pyvista as pv  # noqa: E402

import z3st  # noqa: E402,F401  (installs the house plot style)
from z3st.utils.non_regression import load_yaml  # noqa: E402
from z3st.utils.plotstyle import OKABE_ITO  # noqa: E402

pv.OFF_SCREEN = True
CASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(CASE_DIR, "output")
BLUE, VERMILLION, BLACK = OKABE_ITO[0], OKABE_ITO[1], "#000000"

geom = load_yaml(CASE_DIR, "geometry.yaml")
bcs = load_yaml(CASE_DIR, "boundary_conditions.yaml")["mechanical"]["plate"]
a, L = float(geom["a"]), float(geom["Lx"])
P = next(float(b["traction"]) for b in bcs if b["type"] == "Neumann" and b["region"] == "xmax")

grid = pv.read(os.path.join(OUT_DIR, "fields.vtu"))

# Drawn on the linear triangles of mesh.msh, with the fields sampled onto it.
_m = meshio.read(os.path.join(CASE_DIR, "mesh.msh"))
_T = np.vstack([c.data for c in _m.cells if c.type == "triangle"])
plate = pv.UnstructuredGrid(np.hstack([np.full((len(_T), 1), 3), _T]).ravel(),
                            np.full(len(_T), pv.CellType.TRIANGLE), _m.points).sample(grid)
S = np.asarray(plate.point_data["Stress (points)"]).reshape(-1, 3, 3)
plate.point_data["sigma_xx (MPa)"] = S[:, 0, 0] / 1e6
plate.point_data["von Mises (MPa)"] = np.asarray(plate.point_data["VonMises (points)"]) / 1e6


def frame(pl, size):
    """Look straight down on the plate, showing [0, size] x [0, size]."""
    pl.set_background("white")
    pl.view_xy()
    pl.camera.focal_point = (size / 2, size / 2, 0.0)
    pl.camera.position = (size / 2, size / 2, 1.0)
    pl.camera.parallel_projection = True
    pl.camera.parallel_scale = 0.62 * size


EDGES = dict(show_edges=True, edge_color="#1f1f1f", line_width=0.6)

# --.. mesh, whole quarter and close-up --..
pl = pv.Plotter(off_screen=True, shape=(1, 2), window_size=(1600, 800), border=False)
for col, (size, title) in enumerate(((L, "quarter plate"), (3 * a, "close-up, 3 hole radii"))):
    pl.subplot(0, col)
    pl.add_mesh(plate, color="#56B4E9", lighting=False, **EDGES)
    pl.add_text(title, font_size=14, color="black")
    frame(pl, size)
pl.screenshot(os.path.join(OUT_DIR, "mesh.png"), return_img=False)
pl.close()
print("[plots] mesh.png")

# --.. sigma_xx and von Mises around the hole --..
near = plate.clip_box((0, 3 * a, 0, 3 * a, -1, 1), invert=False)
pl = pv.Plotter(off_screen=True, shape=(1, 2), window_size=(1600, 800), border=False)
for col, name in enumerate(("sigma_xx (MPa)", "von Mises (MPa)")):
    pl.subplot(0, col)
    pl.add_mesh(near, scalars=name, cmap="viridis", lighting=False,
                show_edges=True, edge_color="#404040", line_width=0.2,
                scalar_bar_args={"title": name, "color": "black", "vertical": False,
                                 "position_x": 0.15, "position_y": 0.02, "width": 0.7,
                                 "height": 0.06, "title_font_size": 20, "label_font_size": 16,
                                 "fmt": "%.0f"})
    pl.add_text(f"{name.split(' (')[0]}, remote tension {P/1e6:.0f} MPa along x",
                font_size=12, color="black")
    frame(pl, 3 * a)
pl.screenshot(os.path.join(OUT_DIR, "stress_field.png"), return_img=False)
pl.close()
print("[plots] stress_field.png")

# --.. against Kirsch: the section and the hole edge --..
x, y = grid.points[:, 0], grid.points[:, 1]
r = np.hypot(x, y)
Sg = np.asarray(grid.point_data["Stress (points)"]).reshape(-1, 3, 3)

sec = (np.abs(x) < 1e-9) & (r <= 5 * a + 1e-9)
rs = np.linspace(a, 5 * a, 300)
edge = np.abs(r - a) < 1e-9
phi = np.arctan2(y[edge], x[edge])
et = np.c_[-np.sin(phi), np.cos(phi), np.zeros_like(phi)]
stt = np.einsum("ni,nij,nj->n", et, Sg[edge], et)
ph = np.linspace(0, np.pi / 2, 300)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
ax1.plot(rs / a, 1 + 0.5 * (a / rs) ** 2 + 1.5 * (a / rs) ** 4, "-", color=BLACK, lw=2,
         marker="", label="Kirsch")
ax1.plot(r[sec] / a, Sg[sec, 0, 0] / P, "o", color=BLUE, ms=4, linestyle="none", label="Z3ST (nodes)")
ax1.axhline(1.0, color=BLACK, lw=0.8, ls=":")
ax1.set_xlabel(r"distance from the hole centre, $r/a$")
ax1.set_ylabel(r"$\sigma_{xx}/\sigma$ on the section $x = 0$")
ax1.legend(frameon=False)
ax1.grid(True, alpha=0.3)

ax2.plot(np.degrees(ph), 1 - 2 * np.cos(2 * ph), "-", color=BLACK, lw=2, marker="",
         label=r"Kirsch, $1 - 2\cos 2\varphi$")
ax2.plot(np.degrees(phi), stt / P, "s", color=VERMILLION, ms=4, linestyle="none",
         label="Z3ST (nodes)")
ax2.axhline(0.0, color=BLACK, lw=0.8)
ax2.set_xlabel(r"angle from the load direction, $\varphi$ (deg)")
ax2.set_ylabel(r"hoop stress $\sigma_{\theta\theta}/\sigma$ at $r = a$")
ax2.set_xticks(range(0, 91, 15))
ax2.legend(frameon=False)
ax2.grid(True, alpha=0.3)
plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "kirsch.png"), bbox_inches="tight")
plt.close(fig)
print("[plots] kirsch.png")

print("[plots] done.")
