#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: spherical_shell

non-regression script
---------------------
"""

import os
import sys
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt

# Add the z3st/utils directory to sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "..", "utils"))
from utils_extract_vtu import extract_temperature, list_fields, extract_principal_stresses
from utils_verification import pass_fail_check, regression_check
from utils_plot import plot_field_along_r_xyz

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
VTU_FILE = os.path.join(CASE_DIR, "output", "fields.vtu")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")

# Parameters
Ri = 2.000     # inner radius (m)
Ro = 2.500     # outer radius (m)
t = Ro - Ri    # thickness (m)

Pi = 0.000     # internal pressure (Pa)
Po = 0.000     # external pressure (Pa)

nu = 0.280     # Poisson's ratio
E  = 1.8e+11   # Young modulus (Pa)

alpha = 1.4e-5 # (1/K)
k  = 18.5      # (W/mK
q0 = 9.0e+4    # (W/m³)
mu = 27.0      # (1/m)

T0 = 491.0     # (K)
Ti = 494.0     # (K)

TOLERANCE = 6e-1

# ANSI colors
GREEN = "\033[92m"
RED = "\033[91m"
BOLD = "\033[1m"
END = "\033[0m"

# --.. ..- .-.. .-.. --- analytical solution --.. ..- .-.. .-.. ---
def analytic_T(x):
    """Analytical temperature profile for exponential heating with fixed T0 and Ti (in a slab)."""
    Lx = Ro - Ri
    term1 = Ti + (T0 - Ti) * (x / Lx)
    term2 = (q0 / (mu**2 * k)) * ((x / Lx) * (np.exp(-mu * Lx) - 1) - (np.exp(-mu * x) - 1))
    return term1 + term2

def sigma_th(x, T_num):
    """Analytical thermal stress profile, reasonable estimation from temperature profile."""
    # Constraints in 1, 2 or 3 directions:
    # c = 0, 1, 2, respectively
    c = 3.0

    T_mean = 3.0 / (Ro**3 - Ri**3) * np.trapezoid(T_num * x**2, x)
    return alpha * E / (1.0 - c*nu) * (T_mean - T_num)

# --.. ..- .-.. .-.. --- checks --.. ..- .-.. .-.. ---
if not os.path.exists(VTU_FILE):
    raise FileNotFoundError(f"[ERROR] VTU file not found: {VTU_FILE}")
print(f"[INFO] Using VTU file: {VTU_FILE}")

# --.. ..- .-.. .-.. --- list fields --.. ..- .-.. .-.. ---
list_fields(VTU_FILE)

# --. mesh bounds --..
grid = pv.read(VTU_FILE)
xmin, xmax, ymin, ymax, zmin, zmax = grid.bounds
print(f"\n[INFO] Grid bounds:")
print(f"\tx ∈ [{xmin:.4e}, {xmax:.4e}]")
print(f"\ty ∈ [{ymin:.4e}, {ymax:.4e}]")
print(f"\tz ∈ [{zmin:.4e}, {zmax:.4e}]")

# --.. ..- .-.. .-.. --- extract field --.. ..- .-.. .-.. ---
x, y, z, T = extract_temperature(VTU_FILE)

x_ref = np.linspace(0, t, 200)
T_ref = analytic_T(x_ref)
r_ref = x_ref + Ri
slenderness = 0.5 * (Ri + Ro) / t

rr, TT = plot_field_along_r_xyz(
    x, y, z,
    T, f"Temperature (K), R/t = {slenderness:.2f}",
    CASE_DIR,
    color="tab:blue",
    average="round",  # "round", "bins", "kernel", "weighted", False
    decimals=2, n_bins=100,
    r_ref=r_ref, f_ref=T_ref, label_ref="Analytical solution (slab)"
)

r, _, _, sigma1, sigma2, sigma3 = extract_principal_stresses(
    VTU_FILE,
    average="round",  # None, "bins", "weighted", "kernel", "round"
    decimals=2,
    return_coords=True,
)

# --.. ..- .-.. .-.. --- analytical comparison --.. ..- .-.. .-.. ---
plt.figure()

plt.plot(r, sigma1, "r-", lw=2, label=r"$\sigma_1$")
plt.plot(r, sigma2, "g-", lw=2, label=r"$\sigma_2$")
plt.plot(r, sigma3, "b-", lw=2, label=r"$\sigma_3$")

# thermal_stress = sigma_th(rr, TT)
# plt.plot(rr, thermal_stress, color="tab:pink", lw=2, label=r"$\sigma_{th}$")

plt.xlabel("r (m)")
plt.ylabel("Stress (Pa)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(CASE_DIR, "output", "stress_comparison.png"), dpi=200)
plt.close()

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
T_ana = analytic_T(rr)
Tmax_num   = float(np.max(TT))
Tmax_ref   = float(np.max(T_ana))
RelErr_Tmax = abs(Tmax_num - Tmax_ref) / (Tmax_ref if Tmax_ref != 0 else 1.0)

sigma1_mean = 3.0 / (Ro**3 - Ri**3) * np.trapezoid(sigma1 * r**2, r)
sigma2_mean = 3.0 / (Ro**3 - Ri**3) * np.trapezoid(sigma2 * r**2, r)
sigma3_mean = 3.0 / (Ro**3 - Ri**3) * np.trapezoid(sigma3 * r**2, r)

errors = {
    "T_max": {
        "numerical": Tmax_num,
        "reference": Tmax_ref,
        "abs_error": abs(Tmax_num - Tmax_ref),
        "rel_error": RelErr_Tmax
    },
    "sigma1_mean": {
        "numerical": sigma1_mean,
        "reference": 0.0,
        "abs_error": abs(sigma1_mean),
        "rel_error": abs(sigma1_mean) / (E if E != 0 else 1.0)
    },
    "sigma2_mean": {
        "numerical": sigma2_mean,
        "reference": 0.0,
        "abs_error": abs(sigma2_mean),
        "rel_error": abs(sigma2_mean) / (E if E != 0 else 1.0)
    },
    "sigma3_mean": {
        "numerical": sigma3_mean,
        "reference": 0.0,
        "abs_error": abs(sigma3_mean),
        "rel_error": abs(sigma3_mean) / (E if E != 0 else 1.0)
    }
}


# --. pass/fail check --..
all_pass = pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)

# --. regression vs gold --..
regression_check(errors, CASE_DIR)

print("\n[INFO] non-regression completed.\n")
