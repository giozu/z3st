#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: V_creep_verification

Verifies implicit Norton creep (incremental variational principle, tangent by
AD — models/creep_model.py) on a CONSTANT-STRESS uniaxial bar.

The axisymmetric bar carries a constant axial traction sigma at uniform T,
backward Euler is exact.
With A(T) = A0·exp(-Q/RT) and uniaxial Norton flow (eps_cr_zz rate = A·sigma^n,
lateral rates -1/2 of that):

  1. total axial strain    u_z(L)/L = sigma/E + A·sigma^n·t
  2. creep part            u_z(L)/L - sigma/E = A·sigma^n·t   (sharper)
  3. radial strain         u_r(Ro)/Ro = -nu·sigma/E - (1/2)·A·sigma^n·t
     (checks the deviatoric flow direction, not just the magnitude)

One figure: axial strain vs time, Z3ST per-step vs the closed form.
"""

import os
import yaml
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from z3st.utils.utils_verification import *
from z3st.utils.utils_extract_xdmf import *

R_GAS = 8.314462618

CASE_DIR = os.path.dirname(__file__)
OUT = os.path.join(CASE_DIR, "output")
OUT_JSON = os.path.join(OUT, "non-regression.json")

XDMF_FILE = os.path.join(OUT, "fields.xdmf")

# --. case parameters --..
with open(os.path.join(CASE_DIR, "geometry.yaml")) as f:
    geom = yaml.safe_load(f)
Ro, Lz = float(geom["Ro"]), float(geom["Lz"])

with open(os.path.join(CASE_DIR, "input.yaml")) as f:
    inp = yaml.safe_load(f)
t_total = float(max(inp["time"]))
n_steps = int(inp["n_steps"]) - 1            # active increments

with open(os.path.join(CASE_DIR, next(iter(inp["materials"].values())))) as f:
    mat = yaml.safe_load(f)
E, nu = float(mat["E"]), float(mat["nu"])
A0, n, Q = float(mat["creep_A0"]), float(mat["creep_n"]), float(mat["creep_Q"])
T = float(mat["T_initial"])

with open(os.path.join(CASE_DIR, "boundary_conditions.yaml")) as f:
    bcs = yaml.safe_load(f)
sigma = float(next(b["traction"] for b in bcs["mechanical"]["clad"]
                   if b["type"] == "Neumann"))

# --. analytical references --..
A = A0 * np.exp(-Q / (R_GAS * T))
EPS_EL = sigma / E
EPS_CR = A * sigma**n * t_total
EPS_ZZ_REF = EPS_EL + EPS_CR
EPS_RR_REF = -nu * EPS_EL - 0.5 * EPS_CR
TOLERANCE = 1e-3

# --. numerical results --..
x, y, _, u = extract_field_xdmf(XDMF_FILE, field_name="Displacement", step_index=-1)
r, z = np.asarray(x), np.asarray(y)
u = np.asarray(u)

top = np.abs(z - Lz) < 1e-4 * Lz
out = np.abs(r - Ro) < 1e-4 * Ro
eps_zz = float(np.mean(u[top, 1])) / Lz
eps_rr = float(np.mean(u[out, 0])) / Ro
eps_cr_num = eps_zz - EPS_EL

print(f"[INFO] A(T) = {A:.4e} Pa^-{n:.0f} s^-1 at T = {T:.0f} K, "
      f"sigma = {sigma/1e6:.1f} MPa, t = {t_total:.2e} s")
print(f"[INFO] eps_zz  : numerical = {eps_zz:.6e}, analytical = {EPS_ZZ_REF:.6e}")
print(f"[INFO] eps_cr  : numerical = {eps_cr_num:.6e}, analytical = {EPS_CR:.6e}")
print(f"[INFO] eps_rr  : numerical = {eps_rr:.6e}, analytical = {EPS_RR_REF:.6e}")

errors = {
    "axial_strain_total": {
        "numerical": eps_zz,
        "reference": EPS_ZZ_REF,
        "abs_error": float(abs(eps_zz - EPS_ZZ_REF)),
        "rel_error": float(abs(eps_zz - EPS_ZZ_REF) / EPS_ZZ_REF),
    },
    "creep_strain": {
        "numerical": eps_cr_num,
        "reference": EPS_CR,
        "abs_error": float(abs(eps_cr_num - EPS_CR)),
        "rel_error": float(abs(eps_cr_num - EPS_CR) / EPS_CR),
    },
    "radial_strain_deviatoric": {
        "numerical": eps_rr,
        "reference": EPS_RR_REF,
        "abs_error": float(abs(eps_rr - EPS_RR_REF)),
        "rel_error": float(abs(eps_rr - EPS_RR_REF) / abs(EPS_RR_REF)),
    },
}

# --. figure: axial strain accumulation vs the closed form --..
try:
    times = np.linspace(0.0, t_total, n_steps + 1)
    eps_steps = []
    for k in range(n_steps + 1):
        uk = extract_field_xdmf(XDMF_FILE, "Displacement", step_index=k,
                                return_coords=False)
        eps_steps.append(float(np.mean(np.asarray(uk)[top, 1])) / Lz)
    t_line = np.linspace(0.0, t_total, 100)
    plt.figure(figsize=(7, 5))
    plt.plot(t_line / 86400.0, (EPS_EL + A * sigma**n * t_line) * 100, "k-",
             lw=2.5, alpha=0.7, label="Analytical  σ/E + A·σⁿ·t")
    plt.plot(np.array(times) / 86400.0, np.array(eps_steps) * 100, "o", ms=9,
             mfc="none", mec="r", mew=2, label="Z3ST (implicit, AD tangent)")
    plt.xlabel("time (days)")
    plt.ylabel("axial strain (%)")
    plt.title("Norton creep at constant stress: backward Euler is exact\n"
              f"σ = {sigma/1e6:.0f} MPa, T = {T:.0f} K, n = {n:.0f}")
    plt.grid(True, ls=":", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "creep_strain_history.png"), dpi=150)
    print("[INFO] creep_strain_history.png saved")
except Exception as e:
    print(f"[WARNING] history plot skipped: {type(e).__name__}: {e}")

pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] non-regression completed.\n")
