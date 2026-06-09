#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: 1_thin_slab_neumann_2D

non-regression script
---------------------
Steady-state 1D slab (Neumann-Dirichlet).

"""

import os
import re
import yaml
import numpy as np
import matplotlib.pyplot as plt

from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_verification import *

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
VTU_FILE = os.path.join(CASE_DIR, "output", "fields.vtu")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")
MATERIAL_FILE = os.path.join(CASE_DIR, "../../materials/vessel_steel_0.yaml")
GEOMETRY_FILE = os.path.join(CASE_DIR, "geometry.yaml")
BC_FILE = os.path.join(CASE_DIR, "boundary_conditions.yaml")
MESH_GEO_FILE = os.path.join(CASE_DIR, "mesh.geo")

# Geometry, material, boundary conditions:
with open(GEOMETRY_FILE, 'r') as f:
    geom_data = yaml.safe_load(f)

Lx = float(geom_data.get('Lx'))
Ly = float(geom_data.get('Ly'))

with open(MATERIAL_FILE, 'r') as f:
    mat_data = yaml.safe_load(f)

E     = float(mat_data.get('E'))
nu    = float(mat_data.get('nu'))
k     = float(mat_data.get('k'))
alpha = float(mat_data.get('alpha'))
rho   = float(mat_data.get('rho'))
mu    = float(mat_data.get('mu_gamma'))
q0    = float(mat_data.get('gamma_heating'))

with open(BC_FILE, 'r') as f:
    bc_data = yaml.safe_load(f)
thermal_list = bc_data.get('thermal', {}).get('steel', [])

Ti = next((bc['temperature'] for bc in thermal_list if bc.get('type') == 'Dirichlet'), None)
qf = next((bc['flux'] for bc in thermal_list if bc.get('type') == 'Neumann'), None)

To = Ti - qf/k * Lx

with open(MESH_GEO_FILE, 'r') as f:
    content = f.read()
ny = int(re.search(r'ny\s*=\s*(\d+);', content).group(1)) - 1


print(f"[INFO] Geometry loaded: Lx = {Lx} m, Ly = {Ly} m")
print(f"[INFO] Material loaded: E = {E:.2e} Pa, nu = {nu}")
print(f"[INFO] BCs loaded: Ti = {Ti} K, qf = {qf} W/m²")
print(f"[INFO] nFE loaded: ny = {ny}")


y_target, mask_tol = Ly/2, Ly/(2*ny)  # m, m, m (plane selection and tolerance)

TOLERANCE = 3e-3  # - (relative tolerance for non-regression tests)

# --.. ..- .-.. .-.. --- analytic functions  --.. ..- .-.. .-.. ---
def analytic_T(x):
    """Analytical temperature profile (slab)."""
    term1 = Ti + (To - Ti) * (x / Lx)
    term2 = (q0 / (mu**2 * k)) * ((x / Lx) * (np.exp(-mu * Lx) - 1) - (np.exp(-mu * x) - 1))
    return term1 + term2


def sigma_th(x, T_num, c=1.0):
    """Thermal stress profile (linear)."""
    T_mean = np.trapezoid(T_num, x) / (x.max() - x.min())
    return alpha * E / (1.0 - c * nu) * (T_mean - T_num)


# --.. ..- .-.. .-.. --- checks --.. ..- .-.. .-.. ---
list_fields(VTU_FILE)

# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
print(f"[INFO] Target y-plane for extraction: y = {y_target:.4e} m")

# Numerical results
# Temperature
x_T, y_T, _, T_all = extract_field(VTU_FILE, field_name="Temperature")
mask = np.abs(y_T - y_target) < mask_tol
sort_idx = np.argsort(x_T[mask])

x_T = x_T[mask][sort_idx]
T = T_all[mask][sort_idx]

# Stress
x_S, y_S, _, S_all = extract_field(VTU_FILE, field_name="Stress (cells)")
mask = np.abs(y_S - y_target) < mask_tol
sort_idx = np.argsort(x_S[mask])

x_s = x_S[mask][sort_idx]

sigma_xx = S_all[mask, 0][sort_idx]
sigma_yy = S_all[mask, 4][sort_idx]

# Analytical results
T_ref = analytic_T(x_T)
sigma_th_ref = sigma_th(x_s, analytic_T(x_s), c=1.0)
max_sigma_T = np.max(sigma_yy)

# Plot
Pa_to_MPa = 1e-6

plt.figure(figsize=(10, 7))

# Stress
ax1 = plt.gca()
ax1.plot(x_s, sigma_xx * Pa_to_MPa, "bo", label=r"Num. $\sigma_{xx}$", markersize=4, alpha=0.6)
ax1.plot(x_s, sigma_yy * Pa_to_MPa, "ro", label=r"Num. $\sigma_{yy}$", markersize=4, alpha=0.6)
ax1.plot(
    x_s,
    sigma_th_ref * Pa_to_MPa,
    "m--",
    label=r"Approx. $\sigma_{th}$ (ref)",
    linewidth=2.0,
    alpha=0.7,
)

ax1.set_xlabel("x (m)", fontsize=12)
ax1.set_ylabel("Stress (MPa)", fontsize=12)
ax1.grid(True, linestyle="--", alpha=0.7)

# Temperature
ax2 = ax1.twinx()
ax2.plot(x_T, T, "ks", label="Num. Temperature", markersize=3, alpha=0.4)
ax2.plot(x_T, T_ref, "k--", label="Ana. Temperature", linewidth=1.0, alpha=0.8)
ax2.set_ylabel("Temperature (K)", fontsize=12)

# Legend
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines + lines2, labels + labels2, loc="best", frameon=True)

plt.title(
    rf"$T_i$ = {Ti-273.15:.0f}°C, $T_o$ = {To-273.15:.0f}°C,  $Tmax$ = {np.max(T)-273.15:.0f}°C, $\sigma_T$ = {max_sigma_T*1e-6:.1f} MPa",
    pad=15,
    fontsize=14,
)
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

errors = {
    "L2_error_T": {
        "numerical": L2_T,
        "reference": 0.0,
        "abs_error": L2_T,
        "rel_error": RelL2_T,
    },
    "Linf_error_T": {
        "numerical": Linf_T,
        "reference": 0.0,
        "abs_error": Linf_T,
        "rel_error": Linf_T / np.mean(np.abs(T_ref)),
    },
    "T_max": {
        "numerical": Tmax_num,
        "reference": Tmax_ref,
        "abs_error": abs(Tmax_num - Tmax_ref),
        "rel_error": RelErr_Tmax,
    },
    "L2_error_sigma_yy": {
        "numerical": float(err_sigma),
        "reference": 0.0,
        "abs_error": float(err_sigma),
        "rel_error": float(err_sigma),
    },
}

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] non-regression completed.\n")
