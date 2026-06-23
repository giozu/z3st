#!/usr/bin/env python3
"""Shared checks for 2D MA-MOX conductivity cases."""

import glob
import json
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import yaml

from z3st.materials.fuel_profiles import olander_plutonium_factor
from z3st.materials.magni_mox_thermal import k_numpy
from z3st.models.gpr_conductivity import GPRConductivity
from z3st.utils.utils_verification import pass_fail_check, regression_check


CASE_DIR = os.path.dirname(__file__)
OUT = os.path.join(CASE_DIR, "output")
OUT_JSON = os.path.join(OUT, "non-regression.json")
TOLERANCE = 3e-2


def _read_yaml(name):
    with open(os.path.join(CASE_DIR, name), "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def _last_vtu():
    files = sorted(glob.glob(os.path.join(OUT, "fields_*.vtu")))
    single = os.path.join(OUT, "fields.vtu")
    return files[-1] if files else single


def _composition(mat):
    om = mat.get("OM", mat.get("O_M", None))
    x = mat.get("x", 0.0 if om is None else 2.0 - float(om))
    return {
        "Pu": float(mat.get("Pu", mat.get("pu", 0.0))),
        "Am": float(mat.get("Am", mat.get("am", 0.0))),
        "Np": float(mat.get("Np", mat.get("np", 0.0))),
        "x": float(x),
        "p": float(mat.get("p", mat.get("porosity", 0.0))),
        "burnup": float(mat.get("burnup", mat.get("bu", 0.0))),
    }


def _conductivity(mat, burnup=0.0, Pu=None):
    comp = _composition(mat)
    if Pu is not None:
        comp["Pu"] = Pu
    comp["burnup"] = burnup
    k_card = mat["k"]
    if isinstance(k_card, dict) and str(k_card.get("type", "")).lower() in ("gpr", "gaussian_process"):
        model_path = k_card.get("model", k_card.get("path"))
        if not os.path.isabs(model_path):
            model_path = os.path.normpath(os.path.join(CASE_DIR, model_path))
        return GPRConductivity(model_path, **comp, mode=k_card.get("mode", "mean")), "Magni+GPR"
    return lambda T: k_numpy(T, **comp), "Magni"


def _kirchhoff_center(k_func, R, T_s, q_vol):
    target = q_vol * R * R / 4.0
    tmax = T_s + 700.0
    while True:
        T = np.linspace(T_s, tmax, 8000)
        k = np.asarray(k_func(T), dtype=float)
        integ = np.concatenate([[0.0], np.cumsum(0.5 * (k[1:] + k[:-1]) * np.diff(T))])
        if integ[-1] > 1.05 * target:
            return float(np.interp(target, integ, T))
        tmax += 700.0
        if tmax > 4500.0:
            raise RuntimeError("Kirchhoff integration range too small")


def _entry(num, ref):
    ref_abs = max(abs(float(ref)), 1.0)
    return {
        "numerical": float(num),
        "reference": float(ref),
        "abs_error": float(abs(num - ref)),
        "rel_error": float(abs(num - ref) / ref_abs),
    }


geom = _read_yaml("geometry.yaml")
inp = _read_yaml("input.yaml")
mat = _read_yaml(inp["materials"]["fuel"])
R = float(geom["Ro"])
Lz = float(geom["Lz"])
lhr = float(inp["lhr"][-1])
q_vol = lhr / (np.pi * R * R)
T_surface = 650.0

grid = pv.read(_last_vtu())
pts = grid.points
r_all = pts[:, 0]
z_all = pts[:, 1]
T_all = np.asarray(grid.point_data["Temperature"]).ravel()
bu_all = np.asarray(grid.point_data.get("Burnup", np.zeros_like(T_all))).ravel()

z_mid = 0.5 * Lz
line = np.isclose(z_all, z_mid, atol=max(Lz * 1e-6, 1e-9))
if not np.any(line):
    line = np.abs(z_all - z_mid) <= np.min(np.abs(z_all - z_mid)) + 1e-12
order = np.argsort(r_all[line])
r = r_all[line][order]
T = T_all[line][order]
bu = bu_all[line][order]

is_olander = str(mat.get("Pu_profile", "")).lower() == "olander"
Pu_line = np.full_like(T, float(mat.get("Pu", 0.0)), dtype=float)
if is_olander:
    coords = pts[line][order]
    Pu_line = Pu_line * olander_plutonium_factor(coords, bu, mat, model=None)

errors = {}
finite = bool(np.all(np.isfinite(T_all)) and np.all(np.isfinite(bu_all)))
errors["finite_fields"] = {"numerical": 0.0 if finite else 1.0, "reference": 0.0,
                           "abs_error": 0.0 if finite else 1.0, "rel_error": 0.0 if finite else 1.0}

log_path = os.path.join(CASE_DIR, "log_z3st.md")
if os.path.exists(log_path):
    hits = re.findall(r"Integrated fissile power in \S+:\s*([0-9.eE+\-]+)", open(log_path).read())
    if hits:
        errors["integrated_power"] = _entry(float(hits[-1]), lhr * Lz)

if not is_olander:
    k_func, k_label = _conductivity(mat)
    T_center_ref = _kirchhoff_center(k_func, R, T_surface, q_vol)
    errors["center_temperature_kirchhoff"] = _entry(float(np.max(T)), T_center_ref)
else:
    k_func, k_label = _conductivity(mat, burnup=float(np.max(bu)), Pu=float(np.max(Pu_line)))

plt.figure(figsize=(7, 5))
plt.plot(r * 1e3, T, "o-", ms=3.5, lw=1.8, label=f"Tmax={np.max(T):.1f} K")
plt.xlabel("Radius (mm)")
plt.ylabel("Temperature (K)")
plt.title(f"{k_label}: mid-plane radial temperature")
plt.grid(True, ls=":", alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT, "radial_temperature_profile.png"), dpi=150)

plt.figure(figsize=(7, 5))
plt.plot(r * 1e3, Pu_line, "C1-o", ms=3.5, lw=1.8, label="Pu fraction")
plt.xlabel("Radius (mm)")
plt.ylabel("Pu fraction")
plt.grid(True, ls=":", alpha=0.6)
ax2 = plt.gca().twinx()
ax2.plot(r * 1e3, bu, "C2-s", ms=3.0, lw=1.5, label="Burnup")
ax2.set_ylabel("Burnup (MWd/kgHM)")
plt.title("Local composition/burnup state")
plt.tight_layout()
plt.savefig(os.path.join(OUT, "radial_state_profile.png"), dpi=150)

summary = {
    "case": os.path.basename(CASE_DIR),
    "model": k_label,
    "olander": is_olander,
    "T_min_K": float(np.min(T_all)),
    "T_max_K": float(np.max(T_all)),
    "burnup_mean_MWd_kgHM": float(np.mean(bu_all)),
    "burnup_max_MWd_kgHM": float(np.max(bu_all)),
    "Pu_min": float(np.min(Pu_line)),
    "Pu_max": float(np.max(Pu_line)),
}
with open(os.path.join(OUT, "case_summary.json"), "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2)

print(f"[INFO] {summary['case']}: {k_label}, Tmax={summary['T_max_K']:.3f} K, "
      f"Bu_max={summary['burnup_max_MWd_kgHM']:.3f} MWd/kgHM")
pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)
print("\n[INFO] non-regression completed.\n")
