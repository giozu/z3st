#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: verification/fuel/creep_sinh

Verifies the sinh (Nabarro-Herring-type) thermal creep law of
models/creep_model.py, ``creep: nabarro_herring_structure``, with a temperature-dependent Young's
modulus E(T) and specific heat cp(T) on the card, on a constant-stress
uniaxial bar at uniform T (backward Euler is exact). E(T) is checked at
the single T of the bar; cp(T) is only parsed, the thermal analysis is
stationary.

Law (SI): eps_cr rate = A0·exp(-Q/(R·T))·sinh(B·sigma/(R·T)) (1/s), sigma in
Pa, T in K, R = 8.314462618 J/(mol·K); A0, Q, B are the 1.4970 values of
Többe (1975) converted by materials/steel_1515Ti_Tdep.sinh_creep_SI.
With uniaxial flow (lateral rates -1/2 of the axial one):

  1. total axial strain    u_z(L)/L = sigma/E(T) + eps_cr_rate·t
  2. creep part            u_z(L)/L - sigma/E(T) = eps_cr_rate·t
  3. radial strain         u_r(Ro)/Ro = -nu·sigma/E(T) - (1/2)·eps_cr_rate·t
  4. the same creep rate from the thesis units (%/h, cal/mol, MPa), through
     materials/steel_1515Ti_Tdep.creep_rate_1_4970 (checks the conversion)

One figure: axial strain vs time, Z3ST per-step vs the closed form.
"""

import os
import yaml
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from z3st.materials.steel_1515Ti_Tdep import E_Sobolev, creep_rate_1_4970
from z3st.utils.non_regression import finish, metric
from z3st.utils.utils_extract_xdmf import extract_field_xdmf

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
nu = float(mat["nu"])
A0, Q, B = float(mat["creep_A0"]), float(mat["creep_Q"]), float(mat["creep_B"])

with open(os.path.join(CASE_DIR, "boundary_conditions.yaml")) as f:
    bcs = yaml.safe_load(f)
sigma = float(next(b["traction"] for b in bcs["mechanical"]["clad"]
                   if b["type"] == "Neumann"))
T = float(next(b["temperature"] for b in bcs["thermal"]["clad"]
               if b["type"] == "Dirichlet"))

# --. analytical references --..
E = float(E_Sobolev(T))                                      # (Pa)
RATE = A0 * np.exp(-Q / (R_GAS * T)) * np.sinh(B * sigma / (R_GAS * T))   # (1/s)
RATE_THESIS = float(creep_rate_1_4970(T, sigma / 1.0e6))     # (1/s)
EPS_EL = sigma / E
EPS_CR = RATE * t_total
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

print(f"[INFO] E(T) = {E/1e9:.3f} GPa, creep rate = {RATE:.4e} 1/s at T = {T:.0f} K, "
      f"sigma = {sigma/1e6:.1f} MPa, t = {t_total:.2e} s")
print(f"[INFO] creep rate from thesis units = {RATE_THESIS:.4e} 1/s")
print(f"[INFO] eps_zz  : numerical = {eps_zz:.6e}, analytical = {EPS_ZZ_REF:.6e}")
print(f"[INFO] eps_cr  : numerical = {eps_cr_num:.6e}, analytical = {EPS_CR:.6e}")
print(f"[INFO] eps_rr  : numerical = {eps_rr:.6e}, analytical = {EPS_RR_REF:.6e}")

errors = {
    "axial_strain_total": metric(eps_zz, EPS_ZZ_REF),
    "creep_strain": metric(eps_cr_num, EPS_CR),
    "radial_strain_deviatoric": metric(eps_rr, EPS_RR_REF),
    "rate_unit_conversion": metric(RATE, RATE_THESIS),
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
    plt.plot(t_line / 86400.0, (EPS_EL + RATE * t_line) * 100, "k-",
             lw=2.5, alpha=0.7, label="Analytical  σ/E(T) + ε̇(T, σ)·t")
    plt.plot(np.array(times) / 86400.0, np.array(eps_steps) * 100, "o", ms=9,
             mfc="none", mec="r", mew=2, label="Z3ST (implicit, AD tangent)")
    plt.xlabel("time (days)")
    plt.ylabel("axial strain (%)")
    plt.title("Nabarro-Herring (sinh) creep, 1.4970, at constant stress\n"
              f"σ = {sigma/1e6:.0f} MPa, T = {T:.0f} K, E(T) = {E/1e9:.1f} GPa")
    plt.grid(True, ls=":", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "creep_strain_history.png"), dpi=150)
    print("[INFO] creep_strain_history.png saved")
except Exception as e:
    print(f"[WARNING] history plot skipped: {type(e).__name__}: {e}")

finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)
