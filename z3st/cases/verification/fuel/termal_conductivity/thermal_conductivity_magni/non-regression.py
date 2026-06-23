#!/usr/bin/env python3
"""Verification of Magni/GPR thermal conductivity in a heated pellet.

For a solid cylinder with uniform volumetric heating and T(R)=T_s, the
Kirchhoff transform gives

    int[T_s, T(r)] k(T) dT = q''' * (R**2 - r**2) / 4.

This script uses that relation as a semi-analytic reference for nonlinear k(T).
"""

import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import yaml

from z3st.materials.magni_mox_thermal import k_numpy
from z3st.models.gpr_conductivity import GPRConductivity
from z3st.utils.utils_verification import pass_fail_check, regression_check


CASE_DIR = os.path.dirname(__file__)
OUT = os.path.join(CASE_DIR, "output")
OUT_JSON = os.path.join(OUT, "non-regression.json")
TOLERANCE = 3.0e-2

_single = os.path.join(OUT, "fields.vtu")
_steps = sorted(glob.glob(os.path.join(OUT, "fields_*.vtu")))
VTU_FILE = _steps[-1] if _steps else _single


def _mat_value(mat, *keys, default=0.0):
    for key in keys:
        if key in mat:
            return float(mat[key])
    return default


def _composition(mat):
    om = mat.get("OM", mat.get("O_M", None))
    x = _mat_value(mat, "x", default=(0.0 if om is None else 2.0 - float(om)))
    return {
        "Pu": _mat_value(mat, "Pu", "pu", default=0.0),
        "Am": _mat_value(mat, "Am", "am", default=0.0),
        "Np": _mat_value(mat, "Np", "np", default=0.0),
        "x": x,
        "p": _mat_value(mat, "p", "porosity", default=0.0),
        "burnup": _mat_value(mat, "burnup", "bu", default=0.0),
    }


def _conductivity(mat):
    comp = _composition(mat)
    k_card = mat["k"]
    if isinstance(k_card, dict):
        k_type = str(k_card.get("type", "")).lower()
        if k_type in ("magni", "magni_mox"):
            return lambda T: k_numpy(T, **comp), "Magni"
        model_path = k_card.get("model", k_card.get("path"))
        if not os.path.isabs(model_path):
            model_path = os.path.normpath(os.path.join(CASE_DIR, model_path))
        model = GPRConductivity(model_path, **comp, mode=k_card.get("mode", "mean"))
        return model, "Magni+GPR"
    return lambda T: k_numpy(T, **comp), "Magni"


def _kirchhoff_profile(k_func, r, R, T_s, q_vol):
    h_target = q_vol * (R * R - r * r) / 4.0
    Tmax = T_s + 600.0
    while True:
        T_grid = np.linspace(T_s, Tmax, 6000)
        k_grid = np.asarray(k_func(T_grid), dtype=float)
        integ = np.concatenate([[0.0], np.cumsum(0.5 * (k_grid[1:] + k_grid[:-1]) * np.diff(T_grid))])
        if integ[-1] > 1.05 * float(np.max(h_target)):
            break
        Tmax += 600.0
        if Tmax > 4500.0:
            raise RuntimeError("Kirchhoff grid did not cover the requested heat load")
    return np.interp(h_target, integ, T_grid), T_grid, k_grid


geom = yaml.safe_load(open(os.path.join(CASE_DIR, "geometry.yaml")))
inp = yaml.safe_load(open(os.path.join(CASE_DIR, "input.yaml")))
mat = yaml.safe_load(open(os.path.join(CASE_DIR, next(iter(inp["materials"].values())))))
bcs = yaml.safe_load(open(os.path.join(CASE_DIR, "boundary_conditions.yaml")))

R = float(geom["Ro"])
L = float(geom["Lz"])
area = np.pi * R**2
lhr = float(inp["lhr"][0])
q_vol = lhr / area
T_s = float([bc["temperature"] for bc in bcs["thermal"]["fuel"] if bc["type"] == "Dirichlet"][0])

k_func, k_label = _conductivity(mat)

grid = pv.read(VTU_FILE)
r_all = grid.points[:, 0]
z_all = grid.points[:, 1]
T_all = np.asarray(grid.point_data["Temperature"]).ravel()
z_mid = z_all[np.argmin(np.abs(z_all - 0.5 * L))]
mask = np.isclose(z_all, z_mid)
order = np.argsort(r_all[mask])
r = r_all[mask][order]
T_num = T_all[mask][order]

T_ref, T_grid, k_grid = _kirchhoff_profile(k_func, r, R, T_s, q_vol)
rise_ref = T_ref - T_s
rise_num = T_num - T_s
profile_l2 = float(np.sqrt(np.mean((T_num - T_ref) ** 2)))
profile_linf = float(np.max(np.abs(T_num - T_ref)))
rise_scale = max(float(np.mean(np.abs(rise_ref))), 1.0)

T_center_num = float(T_num[0])
T_center_ref = float(T_ref[0])

print(f"[INFO] conductivity model : {k_label}")
print(f"[INFO] q'''               : {q_vol:.6e} W/m3")
print(f"[INFO] T_center numerical : {T_center_num:.3f} K")
print(f"[INFO] T_center reference : {T_center_ref:.3f} K")
print(f"[INFO] radial L2 error    : {profile_l2:.3e} K")

impact = {
    "model": k_label,
    "T_center_reference_K": T_center_ref,
    "T_center_numerical_K": T_center_num,
}
if k_label == "Magni+GPR":
    comp = _composition(mat)
    T_magni, _, _ = _kirchhoff_profile(lambda T: k_numpy(T, **comp), np.array([0.0]), R, T_s, q_vol)
    impact["T_center_magni_reference_K"] = float(T_magni[0])
    impact["T_center_shift_vs_magni_K"] = float(T_center_ref - T_magni[0])
    impact["T_center_shift_vs_magni_percent"] = float(100.0 * (T_center_ref - T_magni[0]) / T_magni[0])
    print(f"[INFO] GPR centre shift vs Magni reference: {impact['T_center_shift_vs_magni_K']:.3f} K")

os.makedirs(OUT, exist_ok=True)
with open(os.path.join(OUT, "conductivity_impact.json"), "w") as f:
    json.dump(impact, f, indent=2)

errors = {
    "center_temperature": {
        "numerical": T_center_num,
        "reference": T_center_ref,
        "abs_error": abs(T_center_num - T_center_ref),
        "rel_error": abs(T_center_num - T_center_ref) / max(abs(T_center_ref - T_s), 1.0),
    },
    "radial_profile_l2": {
        "numerical": profile_l2,
        "reference": 0.0,
        "abs_error": profile_l2,
        "rel_error": profile_l2 / rise_scale,
    },
    "radial_profile_linf": {
        "numerical": profile_linf,
        "reference": 0.0,
        "abs_error": profile_linf,
        "rel_error": profile_linf / max(float(np.max(np.abs(rise_ref))), 1.0),
    },
}

plt.figure(figsize=(7, 5))
plt.plot(r * 1e3, T_ref, "k-", lw=2.0, label="Kirchhoff reference")
plt.plot(r * 1e3, T_num, "ro", ms=4, fillstyle="none", label="Z3ST")
plt.xlabel("radius r (mm)")
plt.ylabel("Temperature (K)")
plt.title(f"{k_label}: radial pellet temperature at z={z_mid*1e3:.1f} mm")
plt.grid(True, ls=":", alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT, "radial_temperature_profile.png"), dpi=160)

plt.figure(figsize=(7, 5))
plt.plot(T_grid, k_numpy(T_grid, **_composition(mat)), "k--", lw=1.8, label="Magni")
plt.plot(T_grid, k_grid, "C1-", lw=2.0, label=k_label)
plt.xlabel("Temperature (K)")
plt.ylabel("k (W/m/K)")
plt.title("Conductivity law sampled by the verification")
plt.grid(True, ls=":", alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT, "conductivity_vs_temperature.png"), dpi=160)

pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)
print("\n[INFO] non-regression completed.\n")
