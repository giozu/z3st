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
import yaml

from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_verification import *

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
VTU_FILE = os.path.join(CASE_DIR, "output", "fields.vtu")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")
OUT_DIR = os.path.join(CASE_DIR, "output")

# Parameters
with open(os.path.join(CASE_DIR, "geometry.yaml")) as f:
    geom = yaml.safe_load(f)
L = float(geom["Lx"])     # cube edge (m)
R = float(geom["ax"])     # cavity radius (m)

with open(os.path.join(CASE_DIR, "input.yaml")) as f:
    inp = yaml.safe_load(f)
mat_path = os.path.join(CASE_DIR, next(iter(inp["materials"].values())))
with open(mat_path) as f:
    mat = yaml.safe_load(f)

Pi = 1.0e6          # internal pressure (Pa)

E = float(mat["E"])       # Young modulus (/)
nu = float(mat["nu"])     # Poisson ratio

TOLERANCE = 2e-1

# --.. ..- .-.. .-.. --- analytical solution --.. ..- .-.. .-.. ---
def sigma_rr_an(r):
    return -Pi * (R / r) ** 3

def sigma_tt_an(r):
    return +0.5 * Pi * (R / r) ** 3

def sigma_vm_an(r):
    return 1.5 * Pi * (R / r) ** 3


# --.. ..- .-.. .-.. --- checks --.. ..- .-.. .-.. ---
list_fields(VTU_FILE)

# --.. ..- .-.. .-.. --- extract field --.. ..- .-.. .-.. ---
x, y, z, S_all = extract_field(VTU_FILE, field_name="Stress (cells)")

r_raw = np.sqrt(x**2 + y**2 + z**2)
theta = np.arccos(np.divide(z, r_raw, out=np.zeros_like(z), where=r_raw > 0))
phi = np.arctan2(y, x)

s_xx, s_yy, s_zz = S_all[:, 0], S_all[:, 4], S_all[:, 8]
s_xy, s_xz, s_yz = S_all[:, 1], S_all[:, 2], S_all[:, 5]

cth, sth = np.cos(theta), np.sin(theta)
cph, sph = np.cos(phi), np.sin(phi)

sigma_rr_raw = (
    s_xx * (sth * cph)**2 +
    s_yy * (sth * sph)**2 +
    s_zz * (cth)**2 +
    2 * s_xy * (sth**2 * sph * cph) +
    2 * s_xz * (sth * cth * cph) +
    2 * s_yz * (sth * cth * sph)
)

sigma_tt_raw = (
    s_xx * (cth * cph)**2 +
    s_yy * (cth * sph)**2 +
    s_zz * (sth)**2 +
    2 * s_xy * (cth**2 * sph * cph) -
    2 * s_xz * (sth * cth * cph) -
    2 * s_yz * (sth * cth * sph)
)

n_bins = 200
r_edges = np.linspace(R, L/2, n_bins + 1)
r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])

sigma_rr_num = np.zeros(n_bins)
sigma_tt_num = np.zeros(n_bins)

for i in range(n_bins):
    mask = (r_raw >= r_edges[i]) & (r_raw < r_edges[i+1])
    if np.any(mask):
        sigma_rr_num[i] = np.mean(sigma_rr_raw[mask])
        sigma_tt_num[i] = np.mean(sigma_tt_raw[mask])

valid = sigma_rr_num != 0
r = r_centers[valid]
sigma_rr_num = sigma_rr_num[valid]
sigma_tt_num = sigma_pp_num = sigma_tt_num[valid]


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

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] non-regression completed.\n")
