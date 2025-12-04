#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: spherical_pressurized_cavity

non-regression script
---------------------

"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

# Add the z3st/utils directory to sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "..", "utils"))
from utils_extract_vtu import extract_spherical_stresses, list_fields
from utils_verification import pass_fail_check, regression_check

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
VTU_FILE = os.path.join(CASE_DIR, "output", "fields.vtu")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")
OUT_DIR = os.path.join(CASE_DIR, "output")
os.makedirs(OUT_DIR, exist_ok=True)

# Parameters
L = 1.0  # cube edge (m)
R = 0.04  # cavity radius (m)
Pi = 1.0e6  # internal pressure (Pa)

E = 2.0e11  # Young modulus (/)
nu = 0.30  # Poisson ratio

TOLERANCE = 2e-1

# ANSI colors
GREEN = "\033[92m"
RED = "\033[91m"
BOLD = "\033[1m"
END = "\033[0m"


# --.. ..- .-.. .-.. --- analytical solution --.. ..- .-.. .-.. ---
def sigma_rr_an(r):
    return -Pi * (R / r) ** 3


def sigma_tt_an(r):
    return +0.5 * Pi * (R / r) ** 3


def sigma_vm_an(r):
    return 1.5 * Pi * (R / r) ** 3


# --.. ..- .-.. .-.. --- checks --.. ..- .-.. .-.. ---
if not os.path.exists(VTU_FILE):
    raise FileNotFoundError(f"[ERROR] VTU file not found at {VTU_FILE}")
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
r_all, _, _, srr_all, stt_all, spp_all = extract_spherical_stresses(
    VTU_FILE, average="bins", decimals=4, n_bins=200
)

# r in [R, L/2]
mask = (r_all >= R) & (r_all <= L / 2.0)
if not np.any(mask):
    raise RuntimeError(f"{RED}[ERROR]{END} No samples found with R <= r <= L/2.")
r = r_all[mask]
sigma_rr_num = srr_all[mask]
sigma_tt_num = stt_all[mask]
sigma_pp_num = spp_all[mask]

sigma_Tresca_num = np.max(
    np.stack(
        [
            np.abs(sigma_rr_num - sigma_tt_num),
            np.abs(sigma_rr_num - sigma_pp_num),
            np.abs(sigma_tt_num - sigma_pp_num),
        ]
    ),
    axis=0,
)

sigma_vm_num = np.sqrt(
    (
        (sigma_rr_num - sigma_tt_num) ** 2
        + (sigma_tt_num - sigma_pp_num) ** 2
        + (sigma_pp_num - sigma_rr_num) ** 2
    )
    / 2.0
)

sigma_rr_ref = sigma_rr_an(r)
sigma_tt_ref = sigma_tt_an(r)
sigma_vm_ref = sigma_vm_an(r)

plt.figure(figsize=(6, 4))

plt.plot(r, sigma_rr_num, "b.", lw=1.5, label=r"$\sigma_{rr}$ (num)")
plt.plot(r, sigma_tt_num, "r.", lw=1.5, label=r"$\sigma_{\theta\theta}$ (num)")
plt.plot(r, sigma_vm_num, "k.", lw=1.5, label=r"$\sigma_\mathrm{vm}$ (num)")
plt.plot(r, sigma_rr_ref, "b--", lw=1.2, label=r"$\sigma_{rr}$ (ana)")
plt.plot(r, sigma_tt_ref, "r--", lw=1.2, label=r"$\sigma_{\theta\theta}$ (ana)")
plt.plot(r, sigma_vm_ref, "k--", lw=1.2, label=r"$\sigma_\mathrm{vm}$ (ana)")

plt.axvline(R, color="blue", lw=0.8, ls=":", label="Cavity radius")
plt.axvline(L * 0.5, color="red", lw=0.8, ls="--", label="Box edge")

plt.xlabel("r (m)")
plt.ylabel("Stress (Pa)")
plt.title("Stress comparison", fontsize=10)
plt.grid(True, which="both", ls=":", lw=0.5)
plt.legend(fontsize=8, frameon=False, loc="upper right")
plt.tight_layout(rect=[0, 0, 1, 1])
plt.savefig(
    os.path.join(OUT_DIR, "stress_comparison.png"), dpi=300, bbox_inches="tight", transparent=False
)
plt.close()


# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
def rel_L2(num, ref):
    n = np.linalg.norm(ref)
    return float(np.linalg.norm(num - ref) / (n + 1e-30))


err_tt = rel_L2(sigma_tt_num, sigma_tt_ref)
err_rr = rel_L2(sigma_rr_num, sigma_rr_ref)
err_vm = rel_L2(sigma_vm_num, sigma_vm_ref)

errors = {
    "sigma_rr": {"rel_error": float(err_rr), "numerical": np.mean(sigma_rr_num)},
    "sigma_tt": {"rel_error": float(err_tt), "numerical": np.mean(sigma_tt_num)},
    "sigma_vm": {"rel_error": float(err_vm), "numerical": np.mean(sigma_vm_num)},
}

# --. pass/fail check --..
all_pass = pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)

# --. regression vs gold --..
regression_check(errors, CASE_DIR)

print("\n[INFO] non-regression completed.\n")
