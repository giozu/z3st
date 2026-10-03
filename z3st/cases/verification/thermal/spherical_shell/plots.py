#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.4.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""Figures of the spherical_shell case, written to output/.

  mesh.png                the hexahedral octant mesh
  temperature_field.png   temperature on the mesh
  source_profile.png      gamma heating q(r), with the radial mesh layers
  temperature_profile.png T(r): Z3ST, exact sphere solution, slab approximation,
                          and the error against the exact solution
  stress_profile.png      sigma_rr and sigma_tt against the exact solution
  displacement_profile.png u_r against the exact solution
  mesh_convergence.png    errors against cell size, if convergence_data.txt
                          exists (written by mesh_convergence.py)

Run after the case: python3 plots.py (Allrun does it).
"""

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pyvista as pv  # noqa: E402

import z3st  # noqa: E402,F401  (installs the house plot style)
from z3st.utils.plotstyle import OKABE_ITO  # noqa: E402

pv.OFF_SCREEN = True
CASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(CASE_DIR, "output")
sys.path.insert(0, CASE_DIR)
import exact  # noqa: E402

BLUE, VERMILLION, GREEN, BLACK = OKABE_ITO[0], OKABE_ITO[1], OKABE_ITO[2], "#000000"
P = exact.parameters(CASE_DIR)
grid = pv.read(os.path.join(OUT_DIR, "fields.vtu"))


def save(fig, name):
    path = os.path.join(OUT_DIR, name)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"[plots] {name}")


def render(mesh, name, **kw):
    """Offscreen isometric render of the octant."""
    pl = pv.Plotter(off_screen=True, window_size=(1400, 1200))
    pl.add_mesh(mesh, **kw)
    pl.set_background("white")
    pl.view_isometric()
    pl.camera.azimuth = 180  # look from the outside of the octant
    pl.reset_camera()
    pl.screenshot(os.path.join(OUT_DIR, name), return_img=False)
    pl.close()
    print(f"[plots] {name}")


# --.. mesh and temperature field --..
# Drawn on the linear hexahedra of mesh.msh: the VTU written by dolfinx stores
# Lagrange cells, which pyvista triangulates when it draws their edges.
import meshio  # noqa: E402

_m = meshio.read(os.path.join(CASE_DIR, "mesh.msh"))
_H = np.vstack([c.data for c in _m.cells if c.type == "hexahedron"])
hexgrid = pv.UnstructuredGrid(np.hstack([np.full((len(_H), 1), 8), _H]).ravel(),
                              np.full(len(_H), pv.CellType.HEXAHEDRON), _m.points)
hexgrid = hexgrid.sample(grid)  # carries Temperature onto the linear mesh
render(hexgrid, "mesh.png", color="#56B4E9", show_edges=True, edge_color="#1f1f1f", line_width=0.8,
       ambient=0.45, diffuse=0.6, specular=0.0)
render(hexgrid, "temperature_field.png", scalars="Temperature", cmap="Oranges", lighting=False,
       show_edges=True, edge_color="#606060", line_width=0.3,
       scalar_bar_args={"title": "T (K)", "color": "black", "vertical": True,
                        "position_x": 0.85, "position_y": 0.2, "height": 0.6,
                        "title_font_size": 22, "label_font_size": 18, "fmt": "%.0f"})

# radii of the nodes and of the cell centres
pts = grid.points
rn = np.linalg.norm(pts, axis=1)
r = np.linspace(P["Ri"], P["Ro"], 600)
layers = np.unique(np.round(rn, 9))

# --.. gamma source and radial mesh layers --..
fig, ax = plt.subplots(figsize=(6.4, 3.6))
ax.plot(r, exact.source(r, P) / 1e6, "-", color=BLUE, lw=2, marker="", label="gamma heating q(r)")
ax.plot(layers, np.zeros_like(layers), "|", color=BLACK, ms=10, mew=1, label="radial mesh layers")
ax.set_xlabel("r (m)")
ax.set_ylabel("q (MW/m$^3$)")
ax.set_xlim(P["Ri"], P["Ro"])
ax.legend(frameon=False)
ax.grid(True, alpha=0.3)
save(fig, "source_profile.png")

# --.. temperature profile and error --..
T = np.asarray(grid.point_data["Temperature"])
Tex_nodes = exact.temperature(rn, P)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.4, 6.2), sharex=True,
                               gridspec_kw={"height_ratios": [2.2, 1]})
ax1.plot(r, exact.temperature(r, P), "-", color=BLACK, lw=2, marker="", label="exact, sphere")
ax1.plot(rn, T, "o", color=BLUE, ms=3, alpha=0.6, mew=0, linestyle="none", label="Z3ST (nodes)", zorder=3)
ax1.plot(r, exact.temperature_slab(r, P), "--", color=VERMILLION, lw=2, marker="",
         label="slab approximation")
ax1.set_ylabel("T (K)")
ax1.legend(frameon=False)
ax1.grid(True, alpha=0.3)
ax2.plot(rn, T - Tex_nodes, "o", color=BLUE, ms=3, alpha=0.6, mew=0, linestyle="none")
ax2.axhline(0.0, color=BLACK, lw=1)
ax2.set_xlabel("r (m)")
ax2.set_ylabel(r"$T - T_\mathrm{exact}$ (K)")
ax2.set_xlim(P["Ri"], P["Ro"])
ax2.grid(True, alpha=0.3)
save(fig, "temperature_profile.png")

# --.. stresses at the cell centres --..
cc = grid.cell_centers().points
rc = np.linalg.norm(cc, axis=1)
nc = cc / rc[:, None]
S = np.asarray(grid.cell_data["Stress (cells)"]).reshape(-1, 3, 3)
srr = np.einsum("ci,cij,cj->c", nc, S, nc)
stt = 0.5 * (np.einsum("cii->c", S) - srr)
_, srr_ex, stt_ex = exact.mechanics(r, P)
fig, ax = plt.subplots(figsize=(6.4, 4.4))
ax.plot(r, stt_ex / 1e6, "-", color=VERMILLION, lw=2, marker="", label=r"$\sigma_{\theta\theta}$ exact")
ax.plot(rc, stt / 1e6, "s", color=VERMILLION, ms=3, alpha=0.6, mew=0, linestyle="none",
        label=r"$\sigma_{\theta\theta}$ Z3ST (cells)")
ax.plot(r, srr_ex / 1e6, "-", color=BLUE, lw=2, marker="", label=r"$\sigma_{rr}$ exact")
ax.plot(rc, srr / 1e6, "o", color=BLUE, ms=3, alpha=0.6, mew=0, linestyle="none",
        label=r"$\sigma_{rr}$ Z3ST (cells)")
ax.axhline(0.0, color=BLACK, lw=0.8)
ax.set_xlabel("r (m)")
ax.set_ylabel("stress (MPa)")
ax.set_xlim(P["Ri"], P["Ro"])
ax.legend(frameon=False, ncol=2, fontsize=10)
ax.grid(True, alpha=0.3)
save(fig, "stress_profile.png")

# --.. radial displacement --..
ur = np.einsum("ni,ni->n", np.asarray(grid.point_data["Displacement"]), pts / rn[:, None])
ur_ex, _, _ = exact.mechanics(r, P)
fig, ax = plt.subplots(figsize=(6.4, 4.0))
ax.plot(rn, ur * 1e3, "o", color=BLUE, ms=3, alpha=0.6, mew=0, linestyle="none", label="Z3ST (nodes)")
ax.plot(r, ur_ex * 1e3, "-", color=BLACK, lw=2, marker="", label="exact")
ax.set_xlabel("r (m)")
ax.set_ylabel(r"$u_r$ (mm)")
ax.set_xlim(P["Ri"], P["Ro"])
ax.legend(frameon=False)
ax.grid(True, alpha=0.3)
save(fig, "displacement_profile.png")

# --.. mesh convergence, if the series has been run --..
data_file = os.path.join(CASE_DIR, "convergence_data.txt")
if os.path.exists(data_file):
    d = np.loadtxt(data_file)
    h = 1.0 / d[:, 1]  # radial cell count sets the scale; every level halves it
    series = [(5, r"$T$ (L2)", BLACK, "o"), (7, r"$\sigma_{rr}$ (L2)", BLUE, "s"),
              (8, r"$\sigma_{\theta\theta}$ (L2)", VERMILLION, "^"), (9, r"$u_r$ (L2)", GREEN, "D")]
    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    for col, lab, c, m in series:
        ax.loglog(h / h[0], d[:, col], "-", color=c, marker=m, ms=7, lw=2, label=lab)
    # slope-2 reference, anchored below the data
    ref = d[0, 9] * 0.5 * (h / h[0]) ** 2
    ax.loglog(h / h[0], ref, ":", color=BLACK, lw=1.5, marker="", label="slope 2")
    ax.set_xlabel(r"relative cell size $h/h_0$")
    ax.set_ylabel("error / reference scale")
    ax.invert_xaxis()
    ax.legend(frameon=False, fontsize=10)
    ax.grid(True, which="both", alpha=0.3)
    save(fig, "mesh_convergence.png")

print("[plots] done.")
