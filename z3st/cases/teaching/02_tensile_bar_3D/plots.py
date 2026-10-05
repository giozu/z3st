#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""Figures of the tensile-bar case, written to output/.

  mesh.png             the quarter bar and its extruded tetrahedral mesh
  stress_field.png     sigma_zz and von Mises on the deformed bar, the undeformed
                       outline behind it; both uniform, both equal to F/A
  strain_field.png     eps_zz, eps_xx and eps_v on the bar deformed x200, in
                       microstrain: stretched along the axis, contracted across it,
                       and a volume that grows
  inclined_plane.png   normal and shear stress on a plane whose normal is at
                       theta from the axis, FE against sigma cos^2 and sigma/2 sin 2
  inclined_strain.png  normal strain and half the engineering shear in a direction
                       at theta from the axis: Mohr's circle of strain

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

bcs = load_yaml(CASE_DIR, "boundary_conditions.yaml")["mechanical"]["steel"]
P = next(float(b["traction"]) for b in bcs if b["type"] == "Neumann" and b["region"] == "top")
WARP = 50  # displacement magnification: the real strain is 1e-3 and invisible

grid = pv.read(os.path.join(OUT_DIR, "fields.vtu"))

# Drawn on the linear tetrahedra of mesh.msh, with the fields sampled onto it.
_m = meshio.read(os.path.join(CASE_DIR, "mesh.msh"))
_T = np.vstack([c.data for c in _m.cells if c.type == "tetra"])
bar = pv.UnstructuredGrid(np.hstack([np.full((len(_T), 1), 4), _T]).ravel(),
                          np.full(len(_T), pv.CellType.TETRA), _m.points).sample(grid)
S = np.asarray(bar.point_data["Stress (points)"]).reshape(-1, 3, 3)
bar.point_data["sigma_zz (MPa)"] = S[:, 2, 2] / 1e6
bar.point_data["von Mises (MPa)"] = np.asarray(bar.point_data["VonMises (points)"]) / 1e6
Eps = np.asarray(bar.point_data["Strain (points)"]).reshape(-1, 3, 3) * 1e6
bar.point_data["eps_zz (1e-6)"] = Eps[:, 2, 2]
bar.point_data["eps_xx (1e-6)"] = Eps[:, 0, 0]
bar.point_data["eps_v (1e-6)"] = np.trace(Eps, axis1=1, axis2=2)


def view(pl):
    pl.set_background("white")
    pl.view_isometric()
    pl.camera.azimuth = 180  # look at the curved side, not at the symmetry planes
    pl.camera.elevation = 10
    pl.reset_camera()
    pl.camera.zoom(1.7)


# --.. mesh --..
pl = pv.Plotter(off_screen=True, window_size=(700, 1300))
pl.add_mesh(bar, color="#56B4E9", show_edges=True, edge_color="#1f1f1f", line_width=0.8,
            ambient=0.45, diffuse=0.6, specular=0.0)
view(pl)
pl.screenshot(os.path.join(OUT_DIR, "mesh.png"), return_img=False)
pl.close()
print("[plots] mesh.png")

# --.. sigma_zz and von Mises on the deformed bar --..
deformed = bar.warp_by_vector("Displacement", factor=WARP)
pl = pv.Plotter(off_screen=True, shape=(1, 2), window_size=(1400, 1300), border=False)
for col, name in enumerate(("sigma_zz (MPa)", "von Mises (MPa)")):
    pl.subplot(0, col)
    pl.add_mesh(bar.extract_feature_edges(), color="#9a9a9a", line_width=1.5)  # undeformed
    pl.add_mesh(deformed, scalars=name, cmap="viridis", clim=(0.0, 250.0),
                show_edges=True, edge_color="#404040", line_width=0.3, lighting=False,
                scalar_bar_args={"title": name, "color": "black", "vertical": True,
                                 "position_x": 0.72, "position_y": 0.2, "height": 0.6,
                                 "title_font_size": 22, "label_font_size": 18, "fmt": "%.0f",
                                 "n_labels": 6})
    pl.add_text(f"{name.split(' (')[0]}, displacement x{WARP}", font_size=14, color="black")
    view(pl)
pl.screenshot(os.path.join(OUT_DIR, "stress_field.png"), return_img=False)
pl.close()
print("[plots] stress_field.png")

# --.. strains on the deformed bar --..
# A diverging map centred on zero: red stretches, blue contracts. x200 rather than x50,
# so that the lateral contraction, 0.3 of the axial strain, is visible too.
deformed = bar.warp_by_vector("Displacement", factor=4 * WARP)
pl = pv.Plotter(off_screen=True, shape=(1, 3), window_size=(1800, 1300), border=False)
for col, name in enumerate(("eps_zz (1e-6)", "eps_xx (1e-6)", "eps_v (1e-6)")):
    pl.subplot(0, col)
    pl.add_mesh(bar.extract_feature_edges(), color="#9a9a9a", line_width=1.5)  # undeformed
    pl.add_mesh(deformed, scalars=name, cmap="RdBu_r", clim=(-1000.0, 1000.0),
                show_edges=True, edge_color="#404040", line_width=0.3, lighting=False,
                scalar_bar_args={"title": name, "color": "black", "vertical": True,
                                 "position_x": 0.72, "position_y": 0.2, "height": 0.6,
                                 "title_font_size": 22, "label_font_size": 18, "fmt": "%.0f",
                                 "n_labels": 5})
    pl.add_text(f"{name.split(' (')[0]}, displacement x{4*WARP}", font_size=14, color="black")
    view(pl)
pl.screenshot(os.path.join(OUT_DIR, "strain_field.png"), return_img=False)
pl.close()
print("[plots] strain_field.png")

# --.. inclined plane: sigma_n and tau against theta --..
sig = np.asarray(grid.cell_data["Stress (cells)"]).reshape(-1, 3, 3).mean(axis=0)
th = np.radians(np.arange(0, 91, 5))
sn, tau = [], []
for a in th:
    n = np.array([np.sin(a), 0.0, np.cos(a)])
    t = sig.T @ n
    sn.append(t @ n)
    tau.append(np.sqrt(max(t @ t - (t @ n) ** 2, 0.0)))
dense = np.radians(np.linspace(0, 90, 400))

fig, ax = plt.subplots(figsize=(6.4, 4.2))
ax.plot(np.degrees(dense), P * np.cos(dense) ** 2 / 1e6, "-", color=BLACK, lw=1.5, marker="",
        label=r"$\sigma\cos^2\theta$")
ax.plot(np.degrees(dense), P * np.sin(2 * dense) / 2e6, "--", color=BLACK, lw=1.5, marker="",
        label=r"$\frac{\sigma}{2}\sin 2\theta$")
ax.plot(np.degrees(th), np.array(sn) / 1e6, "o", color=BLUE, ms=6, linestyle="none",
        label=r"Z3ST $\sigma_n$")
ax.plot(np.degrees(th), np.array(tau) / 1e6, "s", color=VERMILLION, ms=6, linestyle="none",
        label=r"Z3ST $\tau$")
ax.set_xlabel(r"angle of the plane normal from the bar axis, $\theta$ (deg)")
ax.set_ylabel("stress (MPa)")
ax.set_xticks(range(0, 91, 15))
ax.legend(frameon=False)
ax.grid(True, alpha=0.3)
fig.savefig(os.path.join(OUT_DIR, "inclined_plane.png"), bbox_inches="tight")
plt.close(fig)
print("[plots] inclined_plane.png")

# --.. the same for strain: eps_n and gamma/2 in a direction at theta from the axis --..
eps = np.asarray(grid.cell_data["Strain (cells)"]).reshape(-1, 3, 3).mean(axis=0)
en, half_g = [], []
for a in th:
    n = np.array([np.sin(a), 0.0, np.cos(a)])
    m = np.array([np.cos(a), 0.0, -np.sin(a)])   # in the same plane, at 90 deg to n
    en.append(n @ eps @ n)
    half_g.append(abs(m @ eps @ n))
ez, er = eps[2, 2], eps[0, 0]

fig, ax = plt.subplots(figsize=(6.4, 4.2))
ax.plot(np.degrees(dense), 1e6 * (ez * np.cos(dense) ** 2 + er * np.sin(dense) ** 2), "-",
        color=BLACK, lw=1.5, marker="", label=r"$\varepsilon_{zz}\cos^2\theta + \varepsilon_{xx}\sin^2\theta$")
ax.plot(np.degrees(dense), 1e6 * (ez - er) * np.sin(2 * dense) / 2, "--", color=BLACK, lw=1.5,
        marker="", label=r"$\frac{\varepsilon_{zz}-\varepsilon_{xx}}{2}\sin 2\theta$")
ax.plot(np.degrees(th), 1e6 * np.array(en), "o", color=BLUE, ms=6, linestyle="none",
        label=r"Z3ST $\varepsilon_n$")
ax.plot(np.degrees(th), 1e6 * np.array(half_g), "s", color=VERMILLION, ms=6, linestyle="none",
        label=r"Z3ST $\gamma/2$")
ax.axhline(0.0, color=BLACK, lw=0.8)
ax.set_xlabel(r"angle of the direction from the bar axis, $\theta$ (deg)")
ax.set_ylabel(r"strain ($10^{-6}$)")
ax.set_xticks(range(0, 91, 15))
ax.legend(frameon=False, fontsize=9)
ax.grid(True, alpha=0.3)
fig.savefig(os.path.join(OUT_DIR, "inclined_strain.png"), bbox_inches="tight")
plt.close(fig)
print("[plots] inclined_strain.png")

print("[plots] done.")
