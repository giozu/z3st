#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: verification/mechanics/uniaxial_tension_nonlinear
"""

import os

import matplotlib.pyplot as plt
import numpy as np

from z3st.utils.non_regression import case_paths, finish, load_case, metric
from z3st.utils.utils_extract_vtu import extract_field, list_fields

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR, VTU_FILE, OUT_JSON = case_paths(__file__)

# Geometry and material
geom, inp, mat = load_case(CASE_DIR)
Lx, Ly, Lz = float(geom["Lx"]), float(geom["Ly"]), float(geom["Lz"])  # m

E = float(mat["E"])  # Pa
y_target, z_target, mask_tol = Ly / 2, Lz / 2, Ly * 0.1  # Extraction line

# Analytical references
SIGMA_REF = 125e6  # Pa
VON_MISES_REF = 125e6  # Pa
UX_REF = SIGMA_REF * Lx / E  # m

TOLERANCE = 1e-2  # Relative tolerance

# --.. ..- .-.. .-.. --- checks --.. ..- .-.. .-.. ---
list_fields(VTU_FILE)

# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
print(f"[INFO] Target plane: y={y_target:.4e}, z={z_target:.4e}")

# Stress
x_s_all, y_s_all, z_s_all, S_all = extract_field(VTU_FILE, field_name="Stress (cells)")
mask_s = (np.abs(y_s_all - y_target) < mask_tol) & (np.abs(z_s_all - z_target) < mask_tol)
sort_idx_s = np.argsort(x_s_all[mask_s])

x_s_line = x_s_all[mask_s][sort_idx_s]
sigma_xx = S_all[mask_s, 0][sort_idx_s]

_, _, _, vm_all = extract_field(VTU_FILE, field_name="VonMises (cells)")
sigma_vm_line = vm_all[mask_s][sort_idx_s]  # Use same mask/sort as Stress cells

# Displacement
x_u_all, y_u_all, z_u_all, u_all = extract_field(VTU_FILE, field_name="Displacement")
mask_u = (np.abs(y_u_all - y_target) < mask_tol) & (np.abs(z_u_all - z_target) < mask_tol)
sort_idx_u = np.argsort(x_u_all[mask_u])

x_u_line = x_u_all[mask_u][sort_idx_u]
ux_line = u_all[mask_u, 0][sort_idx_u]

# Plot displacement along X
plt.figure(figsize=(10, 5))
plt.plot(x_u_line, ux_line * 1e3, "o-", color="#009E73", label="Z3ST (ux)", markersize=4)
plt.scatter([0, Lx], [0, UX_REF*1e3], color="#D55E00", label="Analytical Ref")
plt.title(f"Displacement u_x profile (y={y_target}, z={z_target})")
plt.xlabel("x-coordinate (m)")
plt.ylabel("u_x (mm)")
plt.grid(True, linestyle="--", alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(CASE_DIR, "output", "plot_ux.png"))

# Plot stresses along X
plt.figure(figsize=(10, 5))
plt.plot(x_s_line, sigma_xx * 1e-6, "s-", color="#E69F00", label="Sigma XX", markersize=4)
plt.plot(x_s_line, sigma_vm_line * 1e-6, "d-", color="#D55E00", label="Von Mises ", markersize=4)
plt.axhline(y=SIGMA_REF * 1e-6, color="k", linestyle=":", label="Analytical Ref")
plt.title(f"Stress profiles (y={y_target}, z={z_target})")
plt.xlabel("x-coordinate (m)")
plt.ylabel("Stress (MPa)")
plt.ylim(124, 126)
plt.grid(True, linestyle="--", alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(CASE_DIR, "output", "plot_stress.png"))

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
sigma_xx_num = np.mean(sigma_xx)
sigma_vm_num = np.mean(sigma_vm_line)
ux_num = np.max(ux_line)

errors = {
    "sigma_xx": metric(sigma_xx_num, SIGMA_REF),
    "sigma_von_mises": metric(sigma_vm_num, VON_MISES_REF),
    "ux_displacement": metric(ux_num, UX_REF),
}

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)
