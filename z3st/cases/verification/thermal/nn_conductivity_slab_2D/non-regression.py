#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: verification/thermal/nn_conductivity_slab_2D

Steady 1D conduction in a slab whose thermal conductivity is a neural network
k = NN(T), trained to fit the law k(T) = 1/(a + b T). For that law the continuum
temperature profile is closed-form (via the Kirchhoff potential), so we verify
the NN-driven solve against the analytic solution.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import yaml

from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_verification import *

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
VTU_FILE = os.path.join(CASE_DIR, "output", "fields.vtu")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")

with open(os.path.join(CASE_DIR, "geometry.yaml")) as f:
    geom = yaml.safe_load(f)
Lx, Ly = float(geom["Lx"]), float(geom["Ly"])  # m

with open(os.path.join(CASE_DIR, "nn_conductor.yaml")) as f:
    mat = yaml.safe_load(f)
a, b = float(mat["k_law"]["a"]), float(mat["k_law"]["b"])  # W/m/K, W/m/K^2

with open(os.path.join(CASE_DIR, "boundary_conditions.yaml")) as f:
    bcs = yaml.safe_load(f)
bc_by_region = {e["region"]: float(e["temperature"]) for e in bcs["thermal"]["slab"]}
T_L, T_R = bc_by_region["xmin"], bc_by_region["xmax"]  # K at x=0 and x=Lx

y_target, mask_tol = Ly / 2, 0.2  # m (extraction line and tolerance)
TOLERANCE = 5e-3  # - (relative tolerance for non-regression tests)


# --.. ..- .-.. .-.. --- analytic solution --.. ..- .-.. .-.. ---
def analytic_T(x):
    # k(T)=1/(a+bT): Kirchhoff potential psi(T) = (1/b) ln(a+bT) is linear in x
    D = np.log(a + b * T_L) / b
    C = (np.log(a + b * T_R) / b - D) / Lx
    return (np.exp(b * (C * x + D)) - a) / b

# --.. ..- .-.. .-.. --- numerical results --.. ..- .-.. .-.. ---
list_fields(VTU_FILE)
print(f"[INFO] Target y-plane for extraction: y = {y_target:.4e} m")

x_T, y_T, z_T, T_all = extract_field(VTU_FILE, field_name="Temperature")
mask = np.abs(y_T - y_target) < mask_tol
sort_idx = np.argsort(x_T[mask])
x_T = x_T[mask][sort_idx]
T = T_all[mask][sort_idx]

T_ref = analytic_T(x_T)

# --.. ..- .-.. .-.. --- plot --.. ..- .-.. .-.. ---
plt.figure(figsize=(8, 5))
xx = np.linspace(0, Lx, 200)
plt.plot(xx, analytic_T(xx), "k-", lw=2, label="analytic  k(T)=1/(a+bT)")
plt.plot(x_T, T, "ro", ms=4, alpha=0.7, label="z3st  k=NN(T) (lagged)")
plt.xlabel("x (m)")
plt.ylabel("Temperature (K)")
plt.title(rf"NN conductivity slab: $T_L$={T_L:.0f} K, $T_R$={T_R:.0f} K")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plot_path = os.path.join(CASE_DIR, "output", "temperature_comparison.png")
plt.savefig(plot_path, dpi=200)
print(f"[INFO] Plot saved in: {plot_path}")

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
L2_T = float(np.sqrt(np.mean((T - T_ref) ** 2)))
Linf_T = float(np.max(np.abs(T - T_ref)))
RelL2_T = float(L2_T / np.mean(np.abs(T_ref)))

Tmax_num, Tmax_ref = float(np.max(T)), float(np.max(T_ref))
RelErr_Tmax = abs(Tmax_num - Tmax_ref) / Tmax_ref

errors = {
    "L2_error_T": {
        "numerical": L2_T, "reference": 0.0,
        "abs_error": L2_T, "rel_error": RelL2_T,
    },
    "Linf_error_T": {
        "numerical": Linf_T, "reference": 0.0,
        "abs_error": Linf_T, "rel_error": Linf_T / np.mean(np.abs(T_ref)),
    },
    "T_max": {
        "numerical": Tmax_num, "reference": Tmax_ref,
        "abs_error": abs(Tmax_num - Tmax_ref), "rel_error": RelErr_Tmax,
    },
}

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] non-regression completed.\n")
