#!/usr/bin/env python3
"""Comparison plots for the Magni and GPR thermal-conductivity verification cases."""

import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import yaml

from z3st.materials.magni_mox_thermal import k_numpy
from z3st.models.gpr_conductivity import GPRConductivity


ROOT = os.path.dirname(__file__)
MAGNI_DIR = os.path.join(ROOT, "thermal_conductivity_magni")
GPR_DIR = os.path.join(ROOT, "thermal_conductivity_GPR")


def _read_yaml(path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def _case_vtu(case_dir):
    output_dir = os.path.join(case_dir, "output")
    steps = sorted(glob.glob(os.path.join(output_dir, "fields_*.vtu")))
    single = os.path.join(output_dir, "fields.vtu")
    return steps[-1] if steps else single


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


def _gpr_model(mat, case_dir):
    k_card = mat["k"]
    model_path = k_card.get("model", k_card.get("path"))
    if not os.path.isabs(model_path):
        model_path = os.path.normpath(os.path.join(case_dir, model_path))
    return GPRConductivity(model_path, **_composition(mat), mode=k_card.get("mode", "mean"))


def _radial_profile(case_dir):
    grid = pv.read(_case_vtu(case_dir))
    pts = grid.points
    r = pts[:, 0]
    z = pts[:, 1]
    T = np.asarray(grid.point_data["Temperature"]).ravel()
    geom = _read_yaml(os.path.join(case_dir, "geometry.yaml"))
    z_mid = 0.5 * float(geom["Lz"])
    line = np.isclose(z, z_mid, atol=max(float(geom["Lz"]) * 1e-6, 1e-9))
    if not np.any(line):
        line = np.abs(z - z_mid) <= np.min(np.abs(z - z_mid)) + 1e-12
    order = np.argsort(r[line])
    return r[line][order], T[line][order]


def main():
    magni_mat = _read_yaml(os.path.join(MAGNI_DIR, "fuel.yaml"))
    gpr_mat = _read_yaml(os.path.join(GPR_DIR, "fuel.yaml"))
    comp = _composition(magni_mat)
    gpr = _gpr_model(gpr_mat, GPR_DIR)

    r_m, T_m = _radial_profile(MAGNI_DIR)
    r_g, T_g = _radial_profile(GPR_DIR)
    t_min = min(float(np.min(T_m)), float(np.min(T_g)))
    t_max = max(float(np.max(T_m)), float(np.max(T_g)))
    T_probe = np.linspace(t_min, t_max, 300)

    plt.figure(figsize=(7, 5))
    plt.plot(T_probe, k_numpy(T_probe, **comp), "k-", lw=2.4, label="Magni")
    plt.plot(T_probe, gpr(T_probe), "C3--", lw=2.4, label="Magni + GPR")
    plt.xlabel("Temperature (K)")
    plt.ylabel("Thermal conductivity (W/m-K)")
    plt.title("Thermal conductivity comparison")
    plt.grid(True, ls=":", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(ROOT, "thermal_conductivity_comparison.png"), dpi=150)

    plt.figure(figsize=(7, 5))
    plt.plot(r_m * 1e3, T_m, "k-o", ms=3.5, lw=1.8, label=f"Magni, Tc={np.max(T_m):.1f} K")
    plt.plot(r_g * 1e3, T_g, "C3-s", ms=3.5, lw=1.8, label=f"Magni + GPR, Tc={np.max(T_g):.1f} K")
    plt.xlabel("Radius (mm)")
    plt.ylabel("Temperature (K)")
    plt.title("Mid-plane radial temperature profile")
    plt.grid(True, ls=":", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(ROOT, "temperature_profile_comparison.png"), dpi=150)

    print("[INFO] thermal_conductivity_comparison.png saved")
    print("[INFO] temperature_profile_comparison.png saved")


if __name__ == "__main__":
    main()
