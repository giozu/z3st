#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: full_cylinder

Full cylinder with uniform volumetric heat generation. Reference is the
analytical radial temperature profile.
"""

import os

import matplotlib.pyplot as plt
import numpy as np

from z3st.utils.non_regression import case_paths, error_metric, finish, load_case, metric
from z3st.utils.utils_extract_xdmf import extract_field_xdmf, list_fields_xdmf

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
XDMF_FILE = os.path.join(CASE_DIR, "output", "fields.xdmf")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")

# Geometry and material
geom, inp, mat = load_case(CASE_DIR)
Ri, Ro, Lz = float(geom["Ri"]), float(geom["Ro"]), float(geom["Lz"])  # m  inner/outer radius, height


Pi, Po = 0.0, 0.0  # Pa         internal and external pressure
k = 2.5  # W/m·K      thermal conductivity
E, nu, alpha = (
    float(mat["E"]),
    float(mat["nu"]),
    float(mat["alpha"]),
)  # Pa, -, 1/K (Young's modulus, Poisson's ratio, thermal expansion)
To = 500.0  # K          outer surface temperature
Lx = Ro - Ri  # m          wall thickness
slenderness = Ri / Lx  # -          slenderness ratio
z_target, z_tol = 0.0, Lz / 20  # m          z-plane for data extraction
LHR = 5.0e2  # linear heat rate (W/m)

TOLERANCE = 2.0e-2  # -          tolerance for non-regression


# --.. ..- .-.. .-.. --- analytical solution --.. ..- .-.. .-.. ---
def analytic_T(r):
    """Analytical temperature profile in a full cylinder with volumetric heat generation."""
    r = np.asarray(r)
    return LHR / (4 * np.pi * k) * (1 - r**2 / Ro**2) + To


def sigma_th(r, T_num, c=1.0):
    """Thermal stress profile (axisymmetric)."""
    T_mean = 2 / (Ro**2 - Ri**2) * np.trapezoid(T_num * r, r)
    return alpha * E / (1.0 - c * nu) * (T_mean - T_num)


def analytical_thermal_stress(x, T):
    """Analytical thermal stress profiles, from pipe equation."""

    sigma_star = alpha * E * (max(T) - min(T)) / (4 * (1.0 - nu))
    sigma_tt = -sigma_star * (1.0 - 3.0 * (x / Ro) ** 2)
    sigma_rr = -sigma_star * (1.0 - (x / Ro) ** 2)
    sigma_zz = sigma_tt + sigma_rr

    return sigma_rr, sigma_tt, sigma_zz


# --.. ..- .-.. .-.. --- checks --.. ..- .-.. .-.. ---
list_fields_xdmf(XDMF_FILE)

# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
print(f"[INFO] Target z-plane for extraction: z = {z_target:.4e} m")

# Numerical results
# Temperature
x_T, z_T, _, T_all = extract_field_xdmf(XDMF_FILE, field_name="Temperature", step_index=-1)
mask = np.abs(z_T - z_target) < z_tol
sort_idx = np.argsort(x_T[mask])

r_T = x_T[mask][sort_idx]
T = T_all[mask][sort_idx]

# Stress
x_S, z_S, _, S_all = extract_field_xdmf(XDMF_FILE, field_name="Stress", step_index=-1)
mask = np.abs(z_S - z_target) < z_tol
sort_idx = np.argsort(x_S[mask])

r_s = x_S[mask][sort_idx]

sigma_rr = S_all[mask, 0][sort_idx]
sigma_tt = S_all[mask, 4][sort_idx]
sigma_zz = S_all[mask, 8][sort_idx]

# Analytical results
sigma_th_ref = sigma_th(r_T, T, c=1.0)
T_ref = analytic_T(r_T)
sigma_rr_ana_th, sigma_tt_ana_th, sigma_zz_ana_th = analytical_thermal_stress(r_s, analytic_T(r_s))

# Numerical maximum thermal stress (hoop)
max_sigma_T = np.max(sigma_tt)

# Plot
Pa_to_MPa = 1e-6

plt.figure(figsize=(10, 7))

# Stress
ax1 = plt.gca()
ax1.plot(
    r_s, sigma_rr * Pa_to_MPa, "o", color="#D55E00", label=r"Num. $\sigma_{rr}$ (Radial)", markersize=4, alpha=0.6
)
ax1.plot(
    r_s, sigma_rr_ana_th * Pa_to_MPa, "-", color="#D55E00", label=r"Ana. $\sigma_{rr}$ (Radial)", linewidth=1.5
)
ax1.plot(
    r_s,
    sigma_tt * Pa_to_MPa,
    "o", color="#009E73",
    label=r"Num. $\sigma_{\theta\theta}$ (Hoop)",
    markersize=4,
    alpha=0.6,
)
ax1.plot(
    r_s,
    sigma_tt_ana_th * Pa_to_MPa,
    "-", color="#009E73",
    label=r"Ana. $\sigma_{\theta\theta}$ (Hoop)",
    linewidth=1.5,
)
ax1.plot(
    r_s, sigma_zz * Pa_to_MPa, "o", color="#0072B2", label=r"Num. $\sigma_{zz}$ (Axial)", markersize=4, alpha=0.6
)
ax1.plot(r_s, sigma_zz_ana_th * Pa_to_MPa, "-", color="#0072B2", label=r"Ana. $\sigma_{zz}$ (Axial)", linewidth=1.5)
ax1.plot(
    r_T,
    sigma_th_ref * Pa_to_MPa,
    "--", color="#CC79A7",
    label=r"Approx. $\sigma_{th}$ (ref)",
    linewidth=2.0,
    alpha=0.7,
)

ax1.set_xlabel("Radius (m)", fontsize=12)
ax1.set_ylabel("Stress (MPa)", fontsize=12)
ax1.grid(True, linestyle="--", alpha=0.7)

# Temperature
ax2 = ax1.twinx()
ax2.plot(r_T, T, "ks", label="Num. Temperature", markersize=3, alpha=0.4)
ax2.plot(r_T, T_ref, "k--", label="Ana. Temperature", linewidth=1.0, alpha=0.8)
ax2.set_ylabel("Temperature (K)", fontsize=12)

# Legend
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines + lines2, labels + labels2, loc="best", frameon=True)

plt.title(
    rf"$Tmax$ = {np.max(T)-273.15:.0f}°C, $R_i/t$ = {slenderness:.2f}, $\sigma_T$ = {max_sigma_T*1e-6:.1f} MPa",
    pad=15,
    fontsize=14,
)
plt.tight_layout()

plot_path = os.path.join(CASE_DIR, "output", "stress_comparison.png")
plt.savefig(plot_path, dpi=300)
print(f"[INFO] Plot saved in: {plot_path}")


# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
err_tt = np.sqrt(np.mean((sigma_tt - sigma_tt_ana_th) ** 2)) / np.sqrt(np.mean(sigma_tt_ana_th**2))
err_rr = np.sqrt(np.mean((sigma_rr - sigma_rr_ana_th) ** 2)) / np.sqrt(np.mean(sigma_rr_ana_th**2))
err_zz = np.sqrt(np.mean((sigma_zz - sigma_zz_ana_th) ** 2)) / np.sqrt(np.mean(sigma_zz_ana_th**2))

L2_T = float(np.sqrt(np.mean((T - T_ref) ** 2)))
Linf_T = float(np.max(np.abs((T - T_ref))))
RelL2_T = float(L2_T / np.mean(np.abs(T_ref)))

Tmax_num = float(np.max(T))
Tmax_ref = float(np.max(T_ref))
RelErr_Tmax = abs(Tmax_num - Tmax_ref) / Tmax_ref

errors = {
    "L2_error_T": error_metric(L2_T, rel=RelL2_T),
    "Linf_error_T": error_metric(Linf_T, rel=Linf_T / np.mean(np.abs(T_ref))),
    "T_max": metric(Tmax_num, Tmax_ref, rel=RelErr_Tmax),
    "L2_error_sigma_rr": error_metric(err_rr),
    "L2_error_sigma_tt": error_metric(err_tt),
    "L2_error_sigma_zz": error_metric(err_zz),
}

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)

print("\n[INFO] Non-regression completed.\n")
