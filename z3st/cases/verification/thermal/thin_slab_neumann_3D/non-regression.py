#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: verification/thermal/thin_slab_neumann_3D

non-regression script
---------------------
Steady-state 3D slab (Neumann-Dirichlet).

"""

import os, re
import yaml
import numpy as np

from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_plot import plotter_sigma_temperature_slab
from z3st.utils.utils_verification import *

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
VTU_FILE = os.path.join(CASE_DIR, "output", "fields.vtu")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")
OUTPUT_DIR = os.path.join(CASE_DIR, "output")
MATERIAL_FILE = os.path.join(CASE_DIR, "../../../../materials/vessel_steel_0.yaml")
GEOMETRY_FILE = os.path.join(CASE_DIR, "geometry.yaml")
MESH_GEO_FILE = os.path.join(CASE_DIR, "mesh.geo")

# Geometry, material, boundary conditions:
with open(GEOMETRY_FILE, 'r') as f:
    geom_data = yaml.safe_load(f)
Lx = float(geom_data.get('Lx'))
Ly = float(geom_data.get('Ly'))
Lz = float(geom_data.get('Lz'))

with open(MATERIAL_FILE, 'r') as f:
    mat_data = yaml.safe_load(f)
E  = float(mat_data.get('E'))
nu = float(mat_data.get('nu'))
k  = float(mat_data.get('k'))
alpha = float(mat_data.get('alpha'))

with open(MESH_GEO_FILE, 'r') as f:
    content = f.read()
nx = float(re.search(rf'nx\s*=\s*([\d\.]+);', content).group(1))
ny = float(re.search(rf'ny\s*=\s*([\d\.]+);', content).group(1))

print(f"[INFO] Geometry loaded: Lx = {Lx} m, Ly = {Ly} m, Lz = {Lz} m")
print(f"[INFO] Material loaded: E = {E:.2e} Pa, nu = {nu}, k = {k:.2e} W/m·K, alpha = {alpha:.2e} 1/K")

Ti, To = 573.0, 583.0  # K (boundary temperature)
y_target, z_target, mask_tol = Ly / 2, Lz / 2, Ly/(ny-1)/2  # m, m, m (plane selection and tolerance)
TOLERANCE = 3e-3  # - (relative tolerance for non-regression tests)

# --.. ..- .-.. .-.. --- analytic functions  --.. ..- .-.. .-.. ---
def analytic_T(x):
    """Analytical temperature profile (slab)."""
    term1 = Ti + (To - Ti) * (x / Lx)
    return term1

def sigma_th(x, T_num, c=1.0):
    """Thermal stress profile (linear)."""
    T_mean = np.trapezoid(T_num, x) / (x.max() - x.min())
    return alpha * E / (1.0 - c * nu) * (T_mean - T_num)

# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
print(f"[INFO] Target y-plane for extraction: y = {y_target:.4e} m")
print(f"[INFO] Target z-plane for extraction: z = {z_target:.4e} m")

# Numerical results
# Temperature
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
sigma_yy = S_all[mask, 4][sort_idx]

# Analytical results
T_ref = analytic_T(x_T)
sigma_th_ref = sigma_th(x_T, T_ref, c=1.0)
max_sigma_T = np.max(sigma_yy)

# Plot
plotter_sigma_temperature_slab(
    x_s=x_s,
    sigma=sigma_yy,
    x_s_ref=x_T,
    sigma_ref=sigma_th_ref,
    T_ref=T_ref,
    x_T=x_T,
    T=T,
    max_sigma_T=max_sigma_T,
    Ti=Ti,
    To=To,
    CASE_DIR=CASE_DIR,
)

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
L2_T = float(np.sqrt(np.mean((T - T_ref) ** 2)))
Linf_T = float(np.max(np.abs((T - T_ref))))
RelL2_T = float(L2_T / np.mean(np.abs(T_ref)))

Tmax_num = float(np.max(T))
Tmax_ref = float(np.max(T_ref))
RelErr_Tmax = abs(Tmax_num - Tmax_ref) / Tmax_ref

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
}

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] non-regression completed.\n")
