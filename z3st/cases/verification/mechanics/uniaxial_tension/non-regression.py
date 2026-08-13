#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: verification/mechanics/uniaxial_tension
"""

import os

import matplotlib.pyplot as plt
import numpy as np

from z3st.utils.non_regression import case_paths, finish, line, load_case, metric

CASE_DIR, VTU_FILE, OUT_JSON = case_paths(__file__)
geom, inp, mat = load_case(CASE_DIR)

Lx, Ly, Lz = float(geom["Lx"]), float(geom["Ly"]), float(geom["Lz"])
E = float(mat["E"])

y_target, z_target, mask_tol = Ly / 2, Lz / 2, Ly * 0.1  # Extraction line

# Analytical references
SIGMA_REF = 125e6  # Pa
VON_MISES_REF = 125e6  # Pa
UX_REF = SIGMA_REF * Lx / E  # m

TOLERANCE = 1e-2  # Relative tolerance

# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
print(f"[INFO] Target plane: y={y_target:.4e}, z={z_target:.4e}")

plane = dict(y=y_target, z=z_target, tol=mask_tol)
x_s_line, sigma_xx = line(VTU_FILE, "Stress (cells)", component=0, **plane)
_, sigma_vm_line = line(VTU_FILE, "VonMises (cells)", **plane)
x_u_line, ux_line = line(VTU_FILE, "Displacement", component=0, **plane)

# Plot displacement along X
plt.figure(figsize=(10, 5))
plt.plot(x_u_line, ux_line * 1e3, "o-", color="#009E73", label=r"z3st u$_{{x}}$", markersize=4)
plt.scatter([0, Lx], [0, UX_REF*1e3], color="#D55E00", label="Analytical")
plt.title(rf"Displacement u$_{{x}}$ profile (y={y_target}, z={z_target})")
plt.xlabel("x-coordinate (m)")
plt.ylabel("u_x (mm)")
plt.grid(True, linestyle="--", alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(CASE_DIR, "output", "plot_ux.png"))

# Plot stresses along X
plt.figure(figsize=(10, 5))
plt.plot(x_s_line, sigma_xx * 1e-6, "s-", color="#E69F00", label=r"$\sigma_{xx}$", markersize=4)
plt.plot(x_s_line, sigma_vm_line * 1e-6, "d-", color="#D55E00", label=r"$\sigma_{VM}$", markersize=4)
plt.axhline(y=SIGMA_REF * 1e-6, color="k", linestyle=":", label="Analytical")
plt.title(rf"Stress profiles (y={y_target}, z={z_target})")
plt.xlabel("x-coordinate (m)")
plt.ylabel("Stress (MPa)")
plt.ylim(124, 126)
plt.grid(True, linestyle="--", alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(CASE_DIR, "output", "plot_stress.png"))

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
errors = {
    "sigma_xx": metric(np.mean(sigma_xx), SIGMA_REF),
    "sigma_von_mises": metric(np.mean(sigma_vm_line), VON_MISES_REF),
    "ux_displacement": metric(np.max(ux_line), UX_REF),
}

finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)
