#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: verification/fuel/axial_table

Verifies the TABULATED axial power source bus (axial_profile ->
fuel_profiles.tabulated_axial -> set_power) on a tall axisymmetric fuel
column. The axial form factor is a piecewise-linear table f(z_i) from the
material card — the standard fuel-performance input (node-wise peaking
factors from a core-physics calculation).

Burnup accumulates as bu(z) = q(z)·t/(rho·HM·8.64e10), so the final burnup
field IS the normalised profile. Three closed-form checks:

  1. accumulation magnitude — nodal-mean burnup = flat closed form
         bu_mean = q_avg * t_total / (rho * HM * 8.64e10)
     (the mean-1 normalisation preserves the average rating);
  2. table-node ratio — bu(z_k)/bu(z_j) = f_k/f_j at two table nodes,
     normalisation-independent and interpolation-exact where mesh nodes
     coincide with table nodes;
  3. peak/mean — max(f) over the trapezoid mean of the table (the continuous
     mean of a piecewise-linear profile is exactly its trapezoid integral).

One figure: the axial burnup profile, Z3ST nodal vs the interpolated table.
"""

import os
import re
import glob
import yaml
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_verification import *

CASE_DIR = os.path.dirname(__file__)
OUT = os.path.join(CASE_DIR, "output")
OUT_JSON = os.path.join(OUT, "non-regression.json")

# Burnup accumulates over the run -> read the *final* step.
_single = os.path.join(OUT, "fields.vtu")
_steps = sorted(glob.glob(os.path.join(OUT, "fields_*.vtu")))
VTU_FILE = _steps[-1] if _steps else _single

# --. geometry + material + history from the case YAML files --..
with open(os.path.join(CASE_DIR, "geometry.yaml")) as f:
    geom = yaml.safe_load(f)
Ro = float(geom["Ro"])
L = float(geom["Lz"])
area = np.pi * Ro**2

with open(os.path.join(CASE_DIR, "input.yaml")) as f:
    inp = yaml.safe_load(f)
lhr = float(inp["lhr"][-1])
t_total = float(max(inp["time"]))

with open(os.path.join(CASE_DIR, next(iter(inp["materials"].values())))) as f:
    mat = yaml.safe_load(f)
rho = float(mat["rho"])
hm = float(mat.get("heavy_metal_fraction", 0.8815))
z_tab = np.asarray(mat["axial_table_z"], dtype=float)
f_tab = np.asarray(mat["axial_table_f"], dtype=float)

# --. analytical references --..
SECONDS_PER_MWD = 8.64e10
q_avg = lhr / area
BU_MEAN_REF = q_avg * t_total / (rho * hm * SECONDS_PER_MWD)

F_MEAN = float(np.trapezoid(f_tab, z_tab) / (z_tab[-1] - z_tab[0]))
PEAK_RATIO_REF = float(f_tab.max() / F_MEAN)
# Table-node ratio: peak node over first node — normalisation cancels.
NODE_RATIO_REF = float(f_tab[2] / f_tab[0])
TOLERANCE = 1e-2

# --. numerical results --..
xr, yz, _, bu = extract_field(VTU_FILE, field_name="Burnup")
z = np.asarray(yz)
bu = np.asarray(bu)

bu_mean = float(np.mean(bu))
bu_at = lambda z0: float(bu[np.argmin(np.abs(z - z0))])
node_ratio = bu_at(z_tab[2]) / max(bu_at(z_tab[0]), 1e-12)
peak_ratio = float(bu.max()) / max(bu_mean, 1e-12)

print(f"[INFO] q_avg = lhr/area = {q_avg:.4e} W/m^3 over t = {t_total:.3e} s")
print(f"[INFO] nodal-mean burnup : numerical = {bu_mean:.6e}, "
      f"analytical = {BU_MEAN_REF:.6e} MWd/kgU")
print(f"[INFO] node ratio f3/f1  : numerical = {node_ratio:.4f}, "
      f"table = {NODE_RATIO_REF:.4f}")
print(f"[INFO] peak/mean ratio   : numerical = {peak_ratio:.4f}, "
      f"analytical (trapezoid) = {PEAK_RATIO_REF:.4f}")

errors = {
    "burnup_mean_closed_form": {
        "numerical": bu_mean,
        "reference": BU_MEAN_REF,
        "abs_error": float(abs(bu_mean - BU_MEAN_REF)),
        "rel_error": float(abs(bu_mean - BU_MEAN_REF) / BU_MEAN_REF),
    },
    "table_node_ratio": {
        "numerical": node_ratio,
        "reference": NODE_RATIO_REF,
        "abs_error": float(abs(node_ratio - NODE_RATIO_REF)),
        "rel_error": float(abs(node_ratio - NODE_RATIO_REF) / NODE_RATIO_REF),
    },
    "peak_mean_ratio": {
        "numerical": peak_ratio,
        "reference": PEAK_RATIO_REF,
        "abs_error": float(abs(peak_ratio - PEAK_RATIO_REF)),
        "rel_error": float(abs(peak_ratio - PEAK_RATIO_REF) / PEAK_RATIO_REF),
    },
}

# --. integrated power (parsed from the solver log) --..
# set_power prints the exact FE integral of the fissile source per step; the
# axial shaping is separable from the 2πr weight, so the integral must equal
# LHR·L (the normalisation preserves the rating).
LOG = os.path.join(CASE_DIR, "log_z3st.md")
P_REF = lhr * L
if os.path.exists(LOG):
    with open(LOG) as f:
        hits = re.findall(r"Integrated fissile power in \S+:\s*([0-9.eE+\-]+)", f.read())
    if hits:
        P_int = float(hits[-1])
        print(f"[INFO] integrated power  : numerical = {P_int:.6e} W, "
              f"analytical LHR·L = {P_REF:.6e} W")
        errors["total_power"] = {
            "numerical": P_int,
            "reference": P_REF,
            "abs_error": float(abs(P_int - P_REF)),
            "rel_error": float(abs(P_int - P_REF) / P_REF),
        }
    else:
        print("[WARNING] no 'Integrated fissile power' line in log — power check skipped.")
else:
    print("[WARNING] log_z3st.md not found — power check skipped.")

# --. figure: axial burnup profile (Z3ST nodal vs interpolated table) --..
try:
    z_line = np.linspace(z_tab[0], z_tab[-1], 200)
    bu_line = BU_MEAN_REF * np.interp(z_line, z_tab, f_tab) / F_MEAN
    plt.figure(figsize=(7, 5))
    plt.plot(bu_line, z_line * 1e2, "k-", lw=2.5, alpha=0.7,
             label="Analytical  bu_mean · f_table(z)/⟨f⟩")
    plt.scatter(bu, z * 1e2, s=10, facecolors="none", edgecolors="r",
                label="Z3ST nodal burnup")
    plt.scatter(BU_MEAN_REF * f_tab / F_MEAN, z_tab * 1e2, marker="s", s=60,
                facecolors="none", edgecolors="b", label="table nodes")
    plt.axvline(BU_MEAN_REF, color="0.5", ls=":", lw=1.2,
                label=f"mean = {BU_MEAN_REF:.1f} MWd/kgU")
    plt.xlabel("burnup (MWd/kgU)")
    plt.ylabel("axial position z (cm)")
    plt.title("Axial burnup distribution: tabulated form factor\n"
              f"peak/mean = {PEAK_RATIO_REF:.3f}")
    plt.grid(True, ls=":", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "burnup_axial_table_profile.png"), dpi=150)
    print("[INFO] burnup_axial_table_profile.png saved")
except Exception as e:
    print(f"[WARNING] axial-profile plot skipped: {type(e).__name__}: {e}")

pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] non-regression completed.\n")
