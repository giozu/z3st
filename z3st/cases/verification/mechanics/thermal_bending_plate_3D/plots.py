#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""Figures of the free-plate case, written to output/.

  mesh.png      the quarter plate and its hexahedral mesh
  deformed.png  the plate deformed (magnified), coloured by u_x, the undeformed
                outline behind it: a dome of curvature kappa = alpha (To - Ti) / Lx
  profiles.png  T, eps_yy and sigma_yy, sigma_zz (per cell) across the thickness at the
                plate centre, FE against the analytic free-plate solution

Run after the case: python3 plots.py (Allrun does it).
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pyvista as pv  # noqa: E402

import z3st  # noqa: E402,F401  (installs the house plot style)
from z3st.utils.non_regression import load_case  # noqa: E402
from z3st.utils.plotstyle import OKABE_ITO  # noqa: E402

pv.OFF_SCREEN = True
CASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(CASE_DIR, "output")
BLUE, VERMILLION, GREEN, BLACK = OKABE_ITO[0], OKABE_ITO[1], OKABE_ITO[2], "#000000"

geom, inp, mat = load_case(CASE_DIR)
Lx = float(geom["Lx"])
E, nu, alpha, T_ref = float(mat["E"]), float(mat["nu"]), float(mat["alpha"]), float(mat["T_ref"])
Ti, To = 600.0, 400.0  # K   as in boundary_conditions.yaml
g = (To - Ti) / Lx
WARP = 10  # displacement magnification

vtu = pv.read(os.path.join(OUT_DIR, "fields.vtu"))
# Order-1 Lagrange hexes are drawn triangulated; as linear hexes the faces show as quads.
grid = pv.UnstructuredGrid(vtu.cells, np.full(vtu.n_cells, pv.CellType.HEXAHEDRON), vtu.points)
grid.point_data.update(vtu.point_data)
grid.point_data["u_x (mm)"] = np.asarray(grid.point_data["Displacement"])[:, 0] * 1e3


def view(pl):
    pl.set_background("white")
    pl.view_isometric()
    pl.reset_camera()
    pl.camera.zoom(1.2)


# --.. mesh --..
pl = pv.Plotter(off_screen=True, window_size=(1100, 900))
pl.add_mesh(grid, color="#56B4E9", show_edges=True, edge_color="#1f1f1f", line_width=0.8,
            ambient=0.45, diffuse=0.6, specular=0.0)
pl.add_axes(color="black")
view(pl)
pl.screenshot(os.path.join(OUT_DIR, "mesh.png"), return_img=False)
pl.close()
print("[plots] mesh.png")

# --.. deformed plate --..
deformed = grid.warp_by_vector("Displacement", factor=WARP)
pl = pv.Plotter(off_screen=True, window_size=(1100, 900))
pl.add_mesh(grid.extract_feature_edges(), color="#9a9a9a", line_width=1.5)  # undeformed
pl.add_mesh(deformed, scalars="u_x (mm)", cmap="viridis", show_edges=True, edge_color="#404040",
            line_width=0.3, lighting=False,
            scalar_bar_args={"title": "u_x (mm)", "color": "black", "vertical": True,
                             "position_x": 0.85, "position_y": 0.2, "height": 0.6,
                             "title_font_size": 22, "label_font_size": 18, "fmt": "%.1f"})
pl.add_text(f"displacement x{WARP}", font_size=14, color="black")
pl.add_axes(color="black")
view(pl)
pl.screenshot(os.path.join(OUT_DIR, "deformed.png"), return_img=False)
pl.close()
print("[plots] deformed.png")

# --.. profiles across the thickness at the plate centre (y = z = 0) --..
pts = grid.points
line = (np.abs(pts[:, 1]) < 1e-9) & (np.abs(pts[:, 2]) < 1e-9)
order = np.argsort(pts[line, 0])
x = pts[line, 0][order]
T = np.asarray(grid.point_data["Temperature"])[line][order]
eps = np.asarray(grid.point_data["Strain (points)"]).reshape(-1, 3, 3)[line][order]

# Stress per cell (the solver's own DG0 value), on the column of cells next to y = z = 0,
# plotted at the cell centres
cc = vtu.cell_centers().points
col = np.isclose(cc[:, 1], cc[:, 1].min(), atol=1e-6) & np.isclose(cc[:, 2], cc[:, 2].min(), atol=1e-6)
order_c = np.argsort(cc[col, 0])
x_c = cc[col, 0][order_c]
sig = np.asarray(vtu.cell_data["Stress (cells)"]).reshape(-1, 3, 3)[col][order_c]

xa = np.linspace(0.0, Lx, 200)
Ta = Ti + g * xa
eps0, kappa = alpha * (0.5 * (Ti + To) - T_ref), alpha * g
mm = 1e3

fig, axs = plt.subplots(1, 3, figsize=(14, 4.2))

ax = axs[0]
ax.plot(xa * mm, Ta - 273.15, "-", color=BLACK, lw=1.5, marker="", label="analytic")
ax.plot(x * mm, T - 273.15, "o", color=VERMILLION, ms=6, linestyle="none", label="Z3ST")
ax.set_ylabel("temperature (°C)")

ax = axs[1]
ax.plot(xa * mm, 1e6 * (eps0 + kappa * (xa - Lx / 2)), "-", color=BLACK, lw=1.5, marker="",
        label=r"$\varepsilon_0 + \kappa\,(x - L/2)$")
ax.plot(xa * mm, 1e6 * alpha * (Ta - T_ref), "--", color=GREEN, lw=1.5, marker="",
        label=r"$\alpha\,\Delta T(x)$")
ax.plot(x * mm, 1e6 * eps[:, 1, 1], "o", color=BLUE, ms=6, linestyle="none", label=r"Z3ST $\varepsilon_{yy}$")
ax.plot(x * mm, 1e6 * eps[:, 2, 2], "x", color=VERMILLION, ms=6, linestyle="none", label=r"Z3ST $\varepsilon_{zz}$")
ax.set_ylabel(r"total strain ($10^{-6}$)")

# Stress-free plate: scale the axis on the stress a held-flat plate would carry,
# so the FE residual is read against something meaningful.
s_flat = alpha * E / (1 - nu) * (0.5 * (Ti + To) - Ta) / 1e6
ax = axs[2]
ax.plot(xa * mm, np.zeros_like(xa), "-", color=BLACK, lw=1.5, marker="", label="analytic, free plate")
ax.plot(xa * mm, s_flat, ":", color=GREEN, lw=1.5, marker="", label="plate held flat (no bending)")
ax.plot(x_c * mm, sig[:, 1, 1] / 1e6, "o", color=BLUE, ms=6, linestyle="none", label=r"Z3ST $\sigma_{yy}$")
ax.plot(x_c * mm, sig[:, 2, 2] / 1e6, "x", color=VERMILLION, ms=6, linestyle="none", label=r"Z3ST $\sigma_{zz}$")
ax.set_ylabel("stress (MPa)")

for ax in axs:
    ax.set_xlabel("x (mm)")
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "profiles.png"), bbox_inches="tight")
plt.close(fig)
print("[plots] profiles.png")

print("[plots] done.")
