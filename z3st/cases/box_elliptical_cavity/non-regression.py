#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: box_elliptical_cavity

non-regression script
---------------------
"""

import os
import sys

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pyvista as pv

# Add the z3st/output directory to sys.path to import extract.py
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "..", "utils"))
from utils_extract_vtu import *
from utils_plot import *
from utils_verification import pass_fail_check, regression_check

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
VTU_FILE = os.path.join(CASE_DIR, "output", "fields.vtu")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")
OUT_DIR = os.path.join(CASE_DIR, "output")

Lx = 0.5  # box length in x (m)
Ly = 0.5  # box length in y (m)
Lz = 0.5  # box length in z (m)

ax = 0.1072  # semi-axis of ellipsoid in x (m)
ay = 0.1072  # semi-axis of ellipsoid in y (m)
az = 0.05  # semi-axis of ellipsoid in z (m)

stress_applied = 1  # Pa

y_target = 0.0  # (m)
z_target = 0.0  # (m)
ex_tol = 5.0e-1

# ANSI colors
GREEN = "\033[92m"
RED = "\033[91m"
BOLD = "\033[1m"
END = "\033[0m"

# --.. ..- .-.. .-.. --- analytical solution --.. ..- .-.. .-.. ---
VON_MISES_REF = 3.5  # (Pa)
SIGMA_XX_REF = 0.7  # (Pa)

TOLERANCE = 6e-2  # relative tolerance for pass/fail

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

# --.. ..- .-.. .-.. --- extract fields --.. ..- .-.. .-.. ---
print(f"[INFO] Target y-plane for extraction: z = {y_target:.4e} m")
print(f"[INFO] Target z-plane for extraction: z = {z_target:.4e} m")

_, _, _, u = extract_displacement(VTU_FILE)
ux, uy, uz = u[:, 0], u[:, 1], u[:, 2]
print(f"[INFO] Displacement extracted: |u|_max = {np.max(np.linalg.norm(u, axis=1)):.3e}")

x_raw, y_raw, z_raw, comps = extract_stress(VTU_FILE, component="all", prefer="cells")
sigma_xx_raw = comps["xx"]

x_vm_raw, y_vm_raw, z_vm_raw, sigma_vm_raw = extract_VonMises(VTU_FILE)

import pandas as pd

decimals = 6

df = pd.DataFrame({"x": x_raw, "sigma_xx": sigma_xx_raw})
df["x_round"] = df["x"].round(decimals)
sigma_xx_avg = df.groupby("x_round", as_index=False)["sigma_xx"].mean()

df_vm = pd.DataFrame({"x": x_vm_raw, "sigma_vm": sigma_vm_raw})
df_vm["x_round"] = df_vm["x"].round(decimals)
sigma_vm_avg = df_vm.groupby("x_round", as_index=False)["sigma_vm"].mean()

df_join = pd.merge(sigma_xx_avg, sigma_vm_avg, on="x_round", how="inner").sort_values("x_round")
x = df_join["x_round"].to_numpy()
y = np.zeros_like(x)
z = np.zeros_like(x)
sigma_xx = df_join["sigma_xx"].to_numpy()
sigma_vm = df_join["sigma_vm"].to_numpy()

# x in [ax, Lx/2]
mask_x = (x >= ax) & (x <= Lx / 2.0)
mask_y = np.abs(y - y_target) < ex_tol
mask_z = np.abs(z - z_target) < ex_tol
mask = mask_x & mask_y & mask_z

if not np.any(mask):
    raise RuntimeError(f"{RED}[ERROR]{END} No samples found with ax <= x <= Lx/2.")

plt.figure(figsize=(6, 4))
plt.plot(x[mask], sigma_xx[mask], "b.", lw=1.5, label=r"$\sigma_{xx}$ (num)")
plt.plot(x[mask], sigma_vm[mask], "r.", lw=1.5, label=r"$\sigma_{vm}$ (num)")

plt.axvline(ax, color="blue", lw=0.8, ls=":", label="Cavity x-semiaxis")
plt.axvline(Lx * 0.5, color="red", lw=0.8, ls="--", label="Box edge")

plt.gca().xaxis.set_major_formatter(ticker.FuncFormatter(lambda X, _: f"{X*1e3:g}"))
plt.gca().yaxis.set_major_formatter(ticker.FuncFormatter(lambda Y, _: f"{Y/stress_applied:g}"))
plt.xlabel("x (mm)")
plt.ylabel(r"Stress / $\sigma_{BCs}$ (/)")
plt.title(r"$\sigma_{BCs}$ = 1 Pa, z-direction", fontsize=10)
plt.grid(True, which="both", ls=":", lw=0.5)
plt.legend(fontsize=8, frameon=False, loc="upper right")
plt.tight_layout(rect=[0, 0, 1, 1])
plt.savefig(
    os.path.join(OUT_DIR, "stress_comparison.png"), dpi=300, bbox_inches="tight", transparent=False
)
plt.close()

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
x_line = x[mask]
idx = np.argsort(x_line)

errors = {
    "sigma_xx": {
        "numerical": float(np.max(sigma_xx[mask][idx])),
        "reference": SIGMA_XX_REF,
        "abs_error": float(abs(np.max(sigma_xx[mask][idx]) - SIGMA_XX_REF)),
        "rel_error": float(abs(np.max(sigma_xx[mask][idx]) - SIGMA_XX_REF) / SIGMA_XX_REF),
    },
    "sigma_von_mises": {
        "numerical": float(np.max(sigma_vm[mask][idx])),
        "reference": VON_MISES_REF,
        "abs_error": float(abs(np.max(sigma_vm[mask][idx]) - VON_MISES_REF)),
        "rel_error": float(abs(np.max(sigma_vm[mask][idx]) - VON_MISES_REF) / VON_MISES_REF),
    },
    "u_displacement": {
        "numerical": float(np.max(np.linalg.norm(u, axis=1))),
        "reference": 0.0,
        "abs_error": 0.0,
        "rel_error": 0.0,
    },
}

# --. pass/fail check --..
all_pass = pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)

# --. regression vs gold --..
regression_check(errors, CASE_DIR)

print("\n[INFO] non-regression completed.\n")
