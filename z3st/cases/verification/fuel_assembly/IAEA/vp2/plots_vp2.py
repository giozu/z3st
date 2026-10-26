#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
IAEA VP2 free thermal bow (non-linear gradient): z3st against the beam
solution and IWGFR/75 Data Table DT1BOW.

    python plots_vp2.py

T is not linear in y, so the curvature is the section average
kappa = alpha * int(T y dA) / int(y^2 dA) (plane sections), here integrated
over the meshed section; u_y(z) = int_0^z kappa(s) (z - s) ds.

Writes output/vp2_profiles.png (T across the section at the top of the core,
deflection along z) and output/vp2_3d.png (deformed duct coloured by T).
"""

import csv
import os

import h5py
import matplotlib.pyplot as plt
import numpy as np
import yaml
from matplotlib.path import Path

os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")
import pyvista as pv  # noqa: E402

from z3st.utils.utils_extract_xdmf import extract_field_xdmf

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "output")
XDMF = os.path.join(OUT, "fields_vp2.xdmf")
REF_CSV = os.path.join(HERE, "reference", "iwgfr75_ex2.csv")

_geom = yaml.safe_load(open(os.path.join(HERE, "geometry.yaml")))
_mat = yaml.safe_load(open(os.path.join(HERE, "steel_iaea.yaml")))
D, H, T_WALL = (float(_geom[k]) for k in ("D", "H", "t"))
ALPHA = float(_mat["alpha"])

# VP2 profile, as in boundary_conditions.yaml
Z_CORE = (1.5, 2.5)                                # core bottom, top (m)
DT_HALF = 25.0                                     # half the hot-cold difference (K)
Y_KNEE = 0.0375                                    # corners 3, 5 (+) and 2, 6 (-) (m)
R_MID = D / 2 - T_WALL / np.sqrt(3)                # mid-wall corner radius, 75 mm (m)


def ramp(z):
    return np.clip((z - Z_CORE[0]) / (Z_CORE[1] - Z_CORE[0]), 0.0, 1.0)


def dT_top(y):
    """Departure from the 525 °C mean at the top of the core (K)."""
    return DT_HALF * np.clip(y / Y_KNEE, -1.0, 1.0)


def kappa_top():
    """alpha int(dT y dA) / int(y^2 dA) over the hexagonal ring of the mesh
    (vertices on +-y, outer corner radius D/2, flats offset by t)."""
    r_o = D / 2
    r_i = r_o - 2 * T_WALL / np.sqrt(3)
    hexagon = lambda r: Path([(r * np.cos(a), r * np.sin(a))
                              for a in np.pi / 2 + np.arange(6) * np.pi / 3])
    g = np.linspace(-r_o, r_o, 3001)
    X, Y = (a.ravel() for a in np.meshgrid(g, g))
    pts = np.column_stack([X, Y])
    ring = hexagon(r_o).contains_points(pts) & ~hexagon(r_i).contains_points(pts)
    y = Y[ring]
    return ALPHA * np.sum(dT_top(y) * y) / np.sum(y**2)


def analytic_uy(z, k_top):
    s = np.linspace(0.0, H, 4001)
    kappa = -k_top * ramp(s)
    return np.array([np.trapezoid(kappa[s <= zi] * (zi - s[s <= zi]), s[s <= zi])
                     for zi in np.atleast_1d(z)])


def references():
    with open(REF_CSV) as f:
        rows = list(csv.DictReader(line for line in f if not line.startswith("#")))
    beam = {float(r["location"].split("z = ")[1].split(" m")[0]): -float(r["value"])
            for r in rows if r["quantity"] == "u_beam"}
    codes = [-float(r["value"]) for r in rows if r["quantity"].startswith("u_codes")]
    return beam, codes


def plot_profiles(x, y, z, U, T, k_top):
    zl = np.unique(np.round(z, 9))
    uy = np.array([U[np.isclose(z, zi), 1].mean() for zi in zl])  # section mean
    zz = np.linspace(0.0, H, 401)
    ua = analytic_uy(zz, k_top)
    iaea, codes = references()

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 5))
    top = np.isclose(z, Z_CORE[1])
    a1.plot(y[top] * 1e3, T[top] - 273.15, ".", color="k", ms=3, label="z3st, nodes at z = 2.5 m")
    yy = np.linspace(-D / 2, D / 2, 301)
    a1.plot(yy * 1e3, 525 + dT_top(yy), "-", color="tab:orange", lw=1.5, label="IWGFR/75 Fig. 1")
    a1.plot(yy * 1e3, 525 + DT_HALF * yy / R_MID, "--", color="0.5", lw=1, label="Example 1 (linear)")
    a1.set_xlabel("y (mm)")
    a1.set_ylabel("T at the top of the core (°C)")
    a1.legend(fontsize=8)
    a1.grid(alpha=0.3)

    a2.plot(uy * 1e3, zl, "o", color="k", ms=3, mfc="none", label="z3st, section mean")
    a2.plot(ua * 1e3, zz, "-", color="tab:orange", lw=1.5, label="beam, κ over the meshed section")
    zi = sorted(iaea)
    a2.plot([iaea[k] for k in zi], zi, "*", color="tab:blue", ms=11, mfc="none",
            label="IWGFR/75 linear theory")
    a2.plot(codes, [H] * 2, "|", color="tab:green", ms=14, mew=2, label="IWGFR/75, range of 11 codes")
    a2.axhspan(*Z_CORE, color="0.9", zorder=0, label="core")
    a2.set_xlabel(r"$u_y$ (mm)")
    a2.set_ylabel("z (m)")
    a2.set_title(f"top: z3st {uy[-1] * 1e3:.2f} mm, IAEA {iaea[H]:.2f} mm "
                 f"({(uy[-1] * 1e3 / iaea[H] - 1) * 100:+.1f} %)", fontsize=10)
    a2.legend(loc="lower left", fontsize=8)
    fig.suptitle(f"IAEA VP2 free thermal bow, non-linear gradient, α = {ALPHA:.2g} 1/K")
    fig.tight_layout()
    path = os.path.join(OUT, "vp2_profiles.png")
    fig.savefig(path, dpi=150)
    print(f"[INFO] {path}")

    for zk in zi:
        u = uy[np.argmin(np.abs(zl - zk))] * 1e3
        print(f"[INFO] z = {zk:.1f} m: z3st {u:.3f} mm, beam {analytic_uy(zk, k_top)[0] * 1e3:.3f} mm, "
              f"IAEA {iaea[zk]:.3f} mm ({(u / iaea[zk] - 1) * 100:+.2f} %)")


def plot_3d(U, T, warp=20.0):
    with h5py.File(XDMF.replace(".xdmf", ".h5"), "r") as f:
        pts = np.array(f["Mesh/mesh/geometry"])
        topo = np.array(f["Mesh/mesh/topology"])
    grid = pv.UnstructuredGrid(
        np.hstack([np.full((len(topo), 1), topo.shape[1]), topo]).ravel(),
        np.full(len(topo), pv.CellType.HEXAHEDRON, dtype=np.uint8), pts)
    grid.point_data["T (°C)"] = T - 273.15
    grid.point_data["u"] = U
    p = pv.Plotter(off_screen=True, shape=(1, 2), window_size=[1400, 1000])
    for col, view in enumerate(("side", "iso")):
        p.subplot(0, col)
        p.add_mesh(grid.extract_feature_edges(), color="gray", line_width=0.5, opacity=0.5)
        p.add_mesh(grid.warp_by_vector("u", factor=warp), scalars="T (°C)", cmap="coolwarm",
                   scalar_bar_args={"title": "T (°C)", "vertical": True, "position_x": 0.85})
        if view == "side":
            p.add_text(f"side view (from +x), deformed x{warp:.0f}\n+y: hot half", font_size=10)
            p.view_yz()
        else:
            p.add_text("isometric", font_size=10)
            p.view_isometric()
        p.add_axes()
        p.reset_camera()
    path = os.path.join(OUT, "vp2_3d.png")
    p.screenshot(path)
    p.close()
    print(f"[INFO] {path}")


if __name__ == "__main__":
    x, y, z, U = extract_field_xdmf(XDMF, "Displacement")
    T = extract_field_xdmf(XDMF, "Temperature", return_coords=False)
    k = kappa_top()
    print(f"[INFO] kappa at the top of the core: {k:.4e} 1/m "
          f"({k / (ALPHA * DT_HALF / R_MID):.4f} x Example 1 on the same section)")
    plot_profiles(x, y, z, U, T, k)
    plot_3d(U, T)
