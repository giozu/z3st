#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
IAEA VP1 free thermal bow: z3st against the beam solution.

    python plots_vp1.py                 # IAEA VP1, this folder
    python plots_vp1.py flats           # VP1 gradient across the flats, flats/
    python plots_vp1.py flats_1515Ti    # 15-15Ti, 50 °C across the flats

With T linear in y across the duct, the curvature is kappa(z) = alpha dT/dy,
and the clamped cantilever deflects as u_y(z) = int_0^z kappa(s) (z - s) ds.
The hot side (+y) expands more, so the duct bows toward -y.

Writes output/vp1_profiles.png (temperatures on the hot and cold side and
centreline deflection along z) and output/vp1_3d.png (deformed duct coloured
by T, side and isometric views) in the case folder.
"""

import csv
import os
import sys

import h5py
import matplotlib.pyplot as plt
import numpy as np
import yaml

os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")
import pyvista as pv  # noqa: E402

from z3st.utils.utils_extract_xdmf import extract_field_xdmf

VARIANT = sys.argv[1] if len(sys.argv) > 1 else ""
HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), VARIANT)
OUT = os.path.join(HERE, "output")

# per variant: field file, material card, half-width HALF_W over which T
# rises by DT_HALF (as in boundary_conditions.yaml), reference tip bow (mm)
_D = 0.150
_A_MID = _D / 2 * np.cos(np.pi / 6) - 0.003 / 2
_CASES = {
    # ANL/NSE-21/30, §2.5 and Table 2.2
    "": ("fields_vp1", "steel_iaea.yaml", _D / 2,
         {"IAEA analytical": -0.4823 * 25.4, "MOOSE (ANL)": -0.4807 * 25.4}),
    "flats": ("fields_vp1_flats", "../steel_iaea.yaml", _D / 2, {}),
    "flats_1515Ti": ("fields_vp1_1515Ti", "15_15Ti_local.yaml", _A_MID, {}),
}
_FIELDS, _MATFILE, HALF_W, REFS = _CASES[VARIANT]
XDMF = os.path.join(OUT, f"{_FIELDS}.xdmf")
REF_CSV = os.path.join(HERE, "reference", "anl_nse-21-30_vp1.csv")  # IAEA VP1 only

_geom = yaml.safe_load(open(os.path.join(HERE, "geometry.yaml")))
_mat = yaml.safe_load(open(os.path.join(HERE, _MATFILE)))

D, H, T_WALL = (float(_geom[k]) for k in ("D", "H", "t"))
A_MID = D / 2 * np.cos(np.pi / 6) - T_WALL / 2    # mid-wall apothem (m)
R_MID = D / 2 - T_WALL / np.sqrt(3)                # mid-wall corner radius (m)
ALPHA = float(_mat["alpha"])

# VP1 profile, as in boundary_conditions.yaml
Z_CORE = (1.5, 2.5)                                # core bottom, top (m)
DT_HALF = 25.0                                     # half the hot-cold difference (K)
Y_MID = A_MID if _geom.get("orientation", "flat") == "flat" else R_MID  # mid-wall on x = 0


def ramp(z):
    return np.clip((z - Z_CORE[0]) / (Z_CORE[1] - Z_CORE[0]), 0.0, 1.0)


def analytic_uy(z):
    """Cantilever deflection under kappa(z) = -alpha * DT_HALF * ramp(z) / HALF_W."""
    s = np.linspace(0.0, H, 4001)
    kappa = -ALPHA * DT_HALF * ramp(s) / HALF_W
    return np.array([np.trapezoid(kappa[s <= zi] * (zi - s[s <= zi]), s[s <= zi])
                     for zi in np.atleast_1d(z)])


def ref_corner_temperatures():
    """ANL/NSE-21/30 Fig. 2.21: corner temperatures (°C) at the top of the core, by corner label."""
    if VARIANT or not os.path.isfile(REF_CSV):
        return {}
    with open(REF_CSV) as f:
        rows = list(csv.DictReader(line for line in f if not line.startswith("#")))
    return {r["location"]: float(r["value"]) for r in rows if r["quantity"] == "T"}


def load():
    x, y, z, U = extract_field_xdmf(XDMF, "Displacement")
    T = extract_field_xdmf(XDMF, "Temperature", return_coords=False)
    return x, y, z, U, T


def plot_profiles(x, y, z, U, T):
    zl = np.unique(np.round(z, 9))
    uy = np.array([U[np.isclose(z, zi), 1].mean() for zi in zl])  # section mean
    zz = np.linspace(0.0, H, 401)
    ua = analytic_uy(zz)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
    for sign, color, label in ((1, "tab:red", "+y (hot)"), (-1, "tab:blue", "-y (cold)")):
        m = (np.abs(x) < 1e-6) & np.isclose(y, sign * Y_MID)
        o = np.argsort(z[m])
        a1.plot(T[m][o] - 273.15, z[m][o], "o", color=color, ms=3, mfc="none",
                label=f"z3st, {label}")
        a1.plot(400 + ramp(zz) * (125 + sign * DT_HALF * Y_MID / HALF_W), zz, "-", color=color, lw=1)
    corners = ref_corner_temperatures()
    for (name, T_c), mk in zip(corners.items(), ("s", "^", "v", "D")):
        a1.plot([T_c] * 2, [Z_CORE[1], H], mk, color="k", ms=6, mfc="none",
                label=f"ANL Fig. 2.21, {name}: {T_c:g} °C")
    if corners:
        a1.plot(400.0, Z_CORE[0], "o", color="k", ms=6, mfc="none", label="ANL, all corners below core: 400 °C")
    a1.axhspan(*Z_CORE, color="0.9", zorder=0, label="core")
    a1.set_xlabel("temperature (°C): z3st mid-wall, ANL outer corners")
    a1.set_ylabel("z (m)")
    a1.legend(loc="upper left", fontsize=8)

    a2.plot(uy * 1e3, zl, "o", color="k", ms=3, mfc="none", label="z3st, section mean")
    a2.plot(ua * 1e3, zz, "-", color="tab:orange", lw=1.5, label="beam, κ = α dT/dy")
    ref_style = list(zip(REFS.items(), ("*", "D"), ("tab:blue", "tab:green")))
    for (name, u), mk, c in ref_style:
        a2.plot(u, H, mk, color=c, ms=9, mfc="none", label=f"{name}, {u:.2f} mm")
    a2.axhspan(*Z_CORE, color="0.9", zorder=0)
    if REFS:  # the tip values differ by < 1 %: zoom on the top of the duct
        ins = a2.inset_axes([0.40, 0.40, 0.56, 0.30])
        ins.plot(uy * 1e3, zl, "o", color="k", ms=4, mfc="none")
        ins.plot(ua * 1e3, zz, "-", color="tab:orange", lw=1.5)
        ins.plot(uy[-1] * 1e3, H, "o", color="k", ms=4)
        tips = [("z3st", uy[-1] * 1e3, "k")] + [(n, u, c) for (n, u), _, c in ref_style]
        for (name, u), mk, c in ref_style:
            ins.plot(u, H, mk, color=c, ms=10, mfc="none")
        for k, (name, u, c) in enumerate(sorted(tips, key=lambda t: t[1])):
            ins.annotate(f"{name.split(' ')[0]}\n{u:.2f}", (u, H), textcoords="offset points",
                         xytext=(0, 10 if k % 2 == 0 else -26), ha="center", fontsize=7, color=c)
        u0 = min(t[1] for t in tips)
        ins.set_xlim(u0 - 0.08, u0 + 0.25)
        ins.set_ylim(H - 0.035, H + 0.02)
        ins.set_title("top of the duct (mm)", fontsize=8)
        ins.tick_params(labelsize=7)
    a2.set_xlabel(r"$u_y$ (mm)")
    a2.set_title(f"top: z3st {uy[-1] * 1e3:.2f} mm, beam {ua[-1] * 1e3:.2f} mm "
                 f"({(uy[-1] / ua[-1] - 1) * 100:+.1f} %)", fontsize=10)
    a2.legend(loc="lower left", fontsize=8)
    fig.suptitle(f"IAEA VP1 free thermal bow{', ' + VARIANT if VARIANT else ''}, α = {ALPHA:.2g} 1/K")
    fig.tight_layout()
    path = os.path.join(OUT, "vp1_profiles.png")
    fig.savefig(path, dpi=150)
    print(f"[INFO] {path}")
    return uy[-1], ua[-1]


def plot_3d(U, T, warp=20.0):
    with h5py.File(XDMF.replace(".xdmf", ".h5"), "r") as f:
        pts = np.array(f["Mesh/mesh/geometry"])
        topo = np.array(f["Mesh/mesh/topology"])
    grid = pv.UnstructuredGrid(
        np.hstack([np.full((len(topo), 1), topo.shape[1]), topo]).ravel(),
        np.full(len(topo), pv.CellType.HEXAHEDRON, dtype=np.uint8), pts)
    grid.point_data["T (°C)"] = T - 273.15
    grid.point_data["u"] = U
    warped = grid.warp_by_vector("u", factor=warp)
    outline = grid.extract_feature_edges()

    p = pv.Plotter(off_screen=True, shape=(1, 2), window_size=[1400, 1000])
    for col, view in enumerate(("side", "iso")):
        p.subplot(0, col)
        p.add_mesh(outline, color="gray", line_width=0.5, opacity=0.5)
        p.add_mesh(warped, scalars="T (°C)", cmap="coolwarm",
                   scalar_bar_args={"title": "T (°C)", "vertical": True,
                                    "position_x": 0.85})
        if view == "side":
            p.add_text(f"side view (from +x), deformed x{warp:.0f}\n"
                       "+y: hot side", font_size=10)
            p.view_yz()
        else:
            p.add_text("isometric", font_size=10)
            p.view_isometric()
        p.add_axes()
        p.reset_camera()
    path = os.path.join(OUT, "vp1_3d.png")
    p.screenshot(path)
    p.close()
    print(f"[INFO] {path}")


if __name__ == "__main__":
    x, y, z, U, T = load()
    fe, beam = plot_profiles(x, y, z, U, T)
    print(f"[INFO] top u_y: z3st {fe * 1e3:.3f} mm, beam {beam * 1e3:.3f} mm")
    plot_3d(U, T)
