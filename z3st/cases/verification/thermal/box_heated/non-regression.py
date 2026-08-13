#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: verification/thermal/box_heated

non-regression script
---------------------
3D box, fully-constrained and heated

"""

import os

import matplotlib.pyplot as plt
import numpy as np

from z3st.utils.non_regression import case_paths, error_metric, finish, load_case, metric
from z3st.utils.utils_extract_vtu import extract_field, list_fields

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR, VTU_FILE, OUT_JSON = case_paths(__file__)

# Geometry and material
geom, inp, mat = load_case(CASE_DIR)
Lx, Ly, Lz = float(geom["Lx"]), float(geom["Ly"]), float(geom["Lz"])  # m (geometry dimensions)

k, E, nu, alpha = (
    float(mat["k"]),
    float(mat["E"]),
    float(mat["nu"]),
    float(mat["alpha"]),
)  # W/m-K, Pa, -, 1/K (thermal conductivity, Young's modulus, Poisson's ratio, thermal expansion)
Ti, To = 500.0, 500.0  # K (boundary temperatures)
y_target, z_target, mask_tol = Ly / 2, Lz / 2, 0.01  # m, m, m (plane selection and tolerance)

TOLERANCE = 1e-3  # - (relative tolerance for non-regression tests)


# --.. ..- .-.. .-.. --- analytic functions  --.. ..- .-.. .-.. ---
def analytic_T(x):
    return Ti * np.ones_like(x)


# --.. ..- .-.. .-.. --- checks --.. ..- .-.. .-.. ---
list_fields(VTU_FILE)

# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
print(f"[INFO] Target y-plane for extraction: y = {y_target:.4e} m")
print(f"[INFO] Target z-plane for extraction: z = {z_target:.4e} m")

# Numerical temperature
x_T, y_T, z_T, T_all = extract_field(VTU_FILE, field_name="Temperature")
mask = (np.abs(y_T - y_target) < mask_tol) & (np.abs(z_T - z_target) < mask_tol)
sort_idx = np.argsort(x_T[mask])

x_T = x_T[mask][sort_idx]
T = T_all[mask][sort_idx]

# Stress
x_S, y_S, z_S, S_all = extract_field(VTU_FILE, field_name="Stress (cells)")
mask = (np.abs(y_S - y_target) < mask_tol) & (np.abs(z_S - z_target) < mask_tol)
sort_idx = np.argsort(x_S[mask])

x_s = x_S[mask][sort_idx]

sigma_xx = S_all[mask, 0][sort_idx]
sigma_xy = S_all[mask, 1][sort_idx]
sigma_xz = S_all[mask, 2][sort_idx]
sigma_yy = S_all[mask, 4][sort_idx]
sigma_zz = S_all[mask, 8][sort_idx]

# Analytical results
T_ref = analytic_T(x_T)
sigma_th_ref = -alpha * E * np.abs(To - 300) / (1 - 2 * nu) * np.ones_like(x_s)
max_sigma_T = -alpha * E * np.abs(To - 300) / (1 - 2 * nu)

# Plot
Pa_to_MPa = 1e-6

plt.figure(figsize=(10, 7))

# Stress
ax1 = plt.gca()
ax1.plot(x_s, sigma_xx * Pa_to_MPa, ".", color="#0072B2", label=r"Num. $\sigma_{xx}$", alpha=0.5)
ax1.plot(x_s, sigma_yy * Pa_to_MPa, ".", color="#009E73", label=r"Num. $\sigma_{yy}$", alpha=0.5)
ax1.plot(x_s, sigma_zz * Pa_to_MPa, ".", color="#56B4E9", label=r"Num. $\sigma_{zz}$", alpha=0.5)
ax1.plot(
    x_s, sigma_th_ref * Pa_to_MPa, "--", color="#CC79A7", label=r"Ana. $\sigma_{th}$", linewidth=2.0, alpha=0.7
)

ax1.set_xlabel("x (m)", fontsize=12)
ax1.set_ylabel("Stress (MPa)", fontsize=12)
ax1.grid(True, linestyle="--", alpha=0.7)
ax1.set_ylim(-1001, -999)

# Temperature
ax2 = ax1.twinx()
ax2.plot(x_T, T, "ks", label="Num. Temperature", markersize=3, alpha=0.4)
ax2.plot(x_T, T_ref, "k--", label="Ana. Temperature", linewidth=1.0, alpha=0.8)
ax2.set_ylabel("Temperature (K)", fontsize=12)
ax2.set_ylim(499, 501)

# Legend
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines + lines2, labels + labels2, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=3)

plt.title(rf"$T_i$ = {Ti:.0f}°C, $\sigma_T$ = {max_sigma_T*1e-6:.1f} MPa", pad=15, fontsize=14)
plt.tight_layout()

plot_path = os.path.join(CASE_DIR, "output", "stress_comparison.png")
plt.savefig(plot_path, dpi=300)
print(f"[INFO] Plot saved in: {plot_path}")

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
L2_T = float(np.sqrt(np.mean((T - T_ref) ** 2)))
Linf_T = float(np.max(np.abs((T - T_ref))))
RelL2_T = float(L2_T / np.mean(np.abs(T_ref)))

Tmax_num = float(np.max(T))
Tmax_ref = float(np.max(T_ref))
RelErr_Tmax = abs(Tmax_num - Tmax_ref) / Tmax_ref

err_sigma = np.sqrt(np.mean((sigma_yy - sigma_th_ref) ** 2)) / np.max(np.abs(sigma_th_ref))
err_shear_xy = np.max(np.abs(sigma_xy)) / np.max(np.abs(sigma_th_ref))
err_shear_xz = np.max(np.abs(sigma_xz)) / np.max(np.abs(sigma_th_ref))

errors = {
    "L2_error_T": error_metric(L2_T, rel=RelL2_T),
    "Linf_error_T": error_metric(Linf_T, rel=Linf_T / np.mean(np.abs(T_ref))),
    "T_max": metric(Tmax_num, Tmax_ref, rel=RelErr_Tmax),
    "L2_error_sigma_yy": error_metric(err_sigma),
    "L2_error_sigma_xy": error_metric(err_shear_xy),
    "L2_error_sigma_xz": error_metric(err_shear_xz),
}

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)
