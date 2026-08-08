#!/usr/bin/env python3
"""Comparison plots for the four 2D MA-MOX conductivity cases."""

import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import yaml

from z3st.materials.fuel_profiles import olander_plutonium_factor
from z3st.materials.magni_mox_thermal import k_numpy
from z3st.models.gpr_conductivity import GPRConductivity


ROOT = os.path.dirname(__file__)
CASES = [
    ("uniform_Pu_porosity_magni", "Uniform Pu,p - Magni", "k", "-"),
    ("uniform_Pu_porosity_GPR", "Uniform Pu,p - GPR", "C3", "--"),
    ("olander_Pu_burnup_magni", "Olander+burnup - Magni", "C0", "-"),
    ("olander_Pu_burnup_GPR", "Olander+burnup - GPR", "C1", "--"),
]


def _read_yaml(path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def _files(case_dir):
    out = os.path.join(case_dir, "output")
    seq = sorted(glob.glob(os.path.join(out, "fields_*.vtu")))
    if seq:
        return seq
    return [os.path.join(out, "fields.vtu")]


def _composition(mat, burnup=0.0, Pu=None):
    om = mat.get("OM", mat.get("O_M", None))
    x = mat.get("x", 0.0 if om is None else 2.0 - float(om))
    return {
        "Pu": float(mat.get("Pu", 0.0) if Pu is None else Pu),
        "Am": float(mat.get("Am", 0.0)),
        "Np": float(mat.get("Np", 0.0)),
        "x": float(x),
        "p": float(mat.get("p", 0.0)),
        "burnup": float(burnup),
    }


def _model(mat, case_dir, burnup=0.0, Pu=None):
    comp = _composition(mat, burnup=burnup, Pu=Pu)
    k_card = mat["k"]
    if isinstance(k_card, dict) and str(k_card.get("type", "")).lower() in ("gpr", "gaussian_process"):
        model_path = k_card.get("model", k_card.get("path"))
        if not os.path.isabs(model_path):
            model_path = os.path.normpath(os.path.join(case_dir, model_path))
        return GPRConductivity(model_path, **comp, mode=k_card.get("mode", "mean"))
    return lambda T: k_numpy(T, **comp)


def _profile(case_dir, vtu):
    geom = _read_yaml(os.path.join(case_dir, "geometry.yaml"))
    mat = _read_yaml(os.path.join(case_dir, "fuel.yaml"))
    grid = pv.read(vtu)
    pts = grid.points
    r = pts[:, 0]
    z = pts[:, 1]
    T = np.asarray(grid.point_data["Temperature"]).ravel()
    bu = np.asarray(grid.point_data.get("Burnup", np.zeros_like(T))).ravel()
    z_mid = 0.5 * float(geom["Lz"])
    line = np.isclose(z, z_mid, atol=max(float(geom["Lz"]) * 1e-6, 1e-9))
    if not np.any(line):
        line = np.abs(z - z_mid) <= np.min(np.abs(z - z_mid)) + 1e-12
    order = np.argsort(r[line])
    Pu = np.full(np.count_nonzero(line), float(mat.get("Pu", 0.0)))
    if str(mat.get("Pu_profile", "")).lower() == "olander":
        Pu = Pu * olander_plutonium_factor(pts[line], bu[line], mat, model=None)
    return r[line][order], T[line][order], bu[line][order], Pu[order]


plt.figure(figsize=(7.2, 5.2))
all_T = []
for case, label, color, ls in CASES:
    case_dir = os.path.join(ROOT, case)
    r, T, bu, Pu = _profile(case_dir, _files(case_dir)[-1])
    all_T.extend(T.tolist())
    plt.plot(r * 1e3, T, color=color, ls=ls, lw=2.0, label=f"{label} ({np.max(T):.0f} K)")
plt.xlabel("Radius (mm)")
plt.ylabel("Temperature (K)")
plt.title("Final mid-plane temperature profiles")
plt.grid(True, ls=":", alpha=0.6)
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(ROOT, "temperature_profile_comparison.png"), dpi=150)

T_grid = np.linspace(max(650.0, min(all_T)), max(all_T), 320)
plt.figure(figsize=(7.2, 5.2))
for case, label, color, ls in CASES:
    case_dir = os.path.join(ROOT, case)
    mat = _read_yaml(os.path.join(case_dir, "fuel.yaml"))
    _, _, bu, Pu = _profile(case_dir, _files(case_dir)[-1])
    model = _model(mat, case_dir, burnup=float(np.max(bu)), Pu=float(np.max(Pu)))
    plt.plot(T_grid, model(T_grid), color=color, ls=ls, lw=2.0, label=label)
plt.xlabel("Temperature (K)")
plt.ylabel("Thermal conductivity (W/m-K)")
plt.title("Effective conductivity sampled at final local state")
plt.grid(True, ls=":", alpha=0.6)
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(ROOT, "thermal_conductivity_comparison.png"), dpi=150)

plt.figure(figsize=(7.2, 5.2))
for case, label, color, ls in CASES:
    case_dir = os.path.join(ROOT, case)
    times = []
    tmax = []
    bumax = []
    for i, vtu in enumerate(_files(case_dir)):
        g = pv.read(vtu)
        times.append(i)
        tmax.append(float(np.max(g.point_data["Temperature"])))
        if "Burnup" in g.point_data:
            bumax.append(float(np.max(g.point_data["Burnup"])))
        else:
            bumax.append(0.0)
    if len(times) > 1:
        plt.plot(bumax, tmax, color=color, ls=ls, marker="o", lw=2.0, label=label)
plt.xlabel("Maximum burnup (MWd/kgHM)")
plt.ylabel("Maximum temperature (K)")
plt.title("Temperature evolution during burnup sweep")
plt.grid(True, ls=":", alpha=0.6)
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(ROOT, "burnup_temperature_evolution.png"), dpi=150)

summary = {}
for case, label, _, _ in CASES:
    path = os.path.join(ROOT, case, "output", "case_summary.json")
    if os.path.exists(path):
        summary[case] = json.load(open(path, encoding="utf-8"))
with open(os.path.join(ROOT, "comparison_summary.json"), "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2)

print("[INFO] temperature_profile_comparison.png saved")
print("[INFO] thermal_conductivity_comparison.png saved")
print("[INFO] burnup_temperature_evolution.png saved")
