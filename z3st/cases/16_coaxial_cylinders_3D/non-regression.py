#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: coaxial_cylinders_3D

non-regression script
---------------------

"""

import os
import yaml
import numpy as np

from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_plot import *
from z3st.utils.utils_verification import *

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
VTU_FILE = os.path.join(CASE_DIR, "output", "fields.vtu")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")
MATERIAL_FILE_1 = os.path.join(CASE_DIR, "../../materials/ceramic.yaml")
MATERIAL_FILE_2 = os.path.join(CASE_DIR, "../../materials/steel.yaml")
GEOMETRY_FILE = os.path.join(CASE_DIR, "geometry.yaml")
MESH_GEO_FILE = os.path.join(CASE_DIR, "mesh.geo")
INPUT_FILE = os.path.join(CASE_DIR, "input.yaml")
BC_FILE = os.path.join(CASE_DIR, "boundary_conditions.yaml")

# Input
with open(INPUT_FILE, 'r') as f:
    input_data = yaml.safe_load(f)

LHR = input_data['lhr'][0]
h_g = input_data['models']['gap_conductance']['value']

# Geometry
with open(GEOMETRY_FILE, 'r') as f:
    geom_data = yaml.safe_load(f)
Ro1 = float(geom_data.get('outer_radius_1'))
Ri2 = float(geom_data.get('inner_radius_2'))
Ro2 = float(geom_data.get('outer_radius_2'))
Lz = float(geom_data.get('Lz'))

# Material
with open(MATERIAL_FILE_1, 'r') as f:
    mat_data_1 = yaml.safe_load(f)

with open(MATERIAL_FILE_2, 'r') as f:
    mat_data_2 = yaml.safe_load(f)

k1 = 2.5 # (W/m-K)
E1 = float(mat_data_1.get('E'))
nu1 = float(mat_data_1.get('nu'))
alpha1 = float(mat_data_1.get('alpha'))

k2 = float(mat_data_2.get('k'))
E2 = float(mat_data_2.get('E'))
nu2 = float(mat_data_2.get('nu'))
alpha2 = float(mat_data_2.get('alpha'))

# Boundary conditions
with open(BC_FILE, 'r') as f:
    bc_data = yaml.safe_load(f)

To2 = bc_data['thermal']['cyl_2'][1]['temperature']

# Target
z_target, z_tol = Lz / 2, Lz / 10  # z-plane for data extraction
TOLERANCE = 5.0e-2  # tolerance for non-regression

# --.. ..- .-.. .-.. --- analytic functions  --.. ..- .-.. .-.. ---
def analytic_T_1(r, T_interface_1):
    r = np.asarray(r)
    return LHR / (4 * np.pi * k1) * (1 - r**2 / Ro1**2) + T_interface_1


def analytic_T_2(r):
    r = np.asarray(r)
    return To2 + LHR / (2 * np.pi * k2) * np.log(Ro2 / r)


def analytic_T_g(r, T_interface_2):
    r = np.asarray(r)
    k_g = h_g * (Ri2 - Ro1)
    return T_interface_2 + LHR / (2 * np.pi * k_g) * np.log(Ri2 / r)


def analytic_T_piecewise(r):
    r = np.asarray(r)

    # interface temperatures
    T_2_at_Ri2 = analytic_T_2(Ri2)
    T_g_at_Ro1 = analytic_T_g(Ro1, T_2_at_Ri2)

    # allocate output
    T_ref = np.zeros_like(r)

    # regions
    mask_1 = r <= Ro1
    mask_g = (r > Ro1) & (r <= Ri2)
    mask_2 = r > Ri2

    T_ref[mask_1] = analytic_T_1(r[mask_1], T_g_at_Ro1)
    T_ref[mask_g] = analytic_T_g(r[mask_g], T_2_at_Ri2)
    T_ref[mask_2] = analytic_T_2(r[mask_2])

    return T_ref

# --.. ..- .-.. .-.. --- checks --.. ..- .-.. .-.. ---
list_fields(VTU_FILE)

# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
print(f"[INFO] Target z-plane for extraction: z = {z_target:.4e} m")

# Numerical results
x_T, y_T, z_T, T_all = extract_temperature(VTU_FILE)
r_T, T = average_section_radial(x_T, y_T, z_T, T_all, z_target=z_target, tol=z_tol, decimals=4)

# r_s, sigma_rr, sigma_tt, sigma_zz = extract_cylindrical_field(
#     filename=VTU_FILE,
#     z_fixed=z_target,
#     tol=z_tol,
#     case_dir=CASE_DIR,
#     field_hint="Stress",
#     data_source="cell",
#     average=True,
#     decimals=6,
# )

# Analytical results
T_ref = analytic_T_piecewise(r_T)

# Plot
plotter_sigma_temperature_cylinder(
    r_s=None,
    r_T=r_T,
    T=T,
    T_ref=T_ref,
    CASE_DIR=CASE_DIR,
)

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
L2_T = float(np.sqrt(np.mean((T - T_ref) ** 2)))
Linf_T = float(np.max(np.abs((T - T_ref))))
RelL2_T = float(L2_T / np.mean(np.abs(T_ref)))

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
}

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] Non-regression completed.\n")
