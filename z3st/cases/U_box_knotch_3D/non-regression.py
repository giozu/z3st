#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: 00_example

non-regression script
---------------------

"""

import os

import numpy as np
import yaml

from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_plot import plot_field_along_x
from z3st.utils.utils_verification import *

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
VTU_FILE = os.path.join(CASE_DIR, "output", "fields.vtu")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")

# Geometry and material
with open(os.path.join(CASE_DIR, "geometry.yaml")) as f:
    geom = yaml.safe_load(f)
Lx, Ly, Lz = float(geom["Lx"]), float(geom["Ly"]), float(geom["Lz"])  # m (geometry dimensions)
with open(os.path.join(CASE_DIR, "input.yaml")) as f:
    inp = yaml.safe_load(f)
mat_path = os.path.join(CASE_DIR, next(iter(inp["materials"].values())))
with open(mat_path) as f:
    mat = yaml.safe_load(f)
E = float(mat["E"])  # (Pa) Young modulus
y_target, z_target, mask_tol = Ly / 2, Lz / 2, 0.01  # m, m, m (plane selection and tolerance)

# --.. ..- .-.. .-.. --- analytic functions  --.. ..- .-.. .-.. ---
SIGMA_REF = 125e6  # (Pa)
VON_MISES_REF = 125e6  # (Pa)
UX_REF = SIGMA_REF * Lx / E  # (m)

TOLERANCE = 1e-2  # relative tolerance for pass/fail

# --.. ..- .-.. .-.. --- checks --.. ..- .-.. .-.. ---
list_fields(VTU_FILE)

# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
print(f"[INFO] Target y-plane for extraction: y = {y_target:.4e} m")
print(f"[INFO] Target z-plane for extraction: z = {z_target:.4e} m")

x_s, y_s, z_s, s = extract_stress(VTU_FILE, component="all", return_coords=True, prefer="cells")
_, sigma_xx = average_section(
    x_s, y_s, z_s, s["xx"], y_target, z_target, mask_tol, decimals=5, label="sigma_xx"
)
_, sigma_yy = average_section(
    x_s, y_s, z_s, s["yy"], y_target, z_target, mask_tol, decimals=5, label="sigma_yy"
)
x_s, sigma_zz = average_section(
    x_s, y_s, z_s, s["zz"], y_target, z_target, mask_tol, decimals=5, label="sigma_zz"
)

x, y, z, u = extract_displacement(VTU_FILE)
ux, uy, uz = u[:, 0], u[:, 1], u[:, 2]
print(f"[INFO] Displacement extracted: |u|_max = {np.max(np.linalg.norm(u, axis=1)):.3e}")

_, _, _, sigma_vm = extract_VonMises(VTU_FILE, prefer="points")

plot_field_along_x(
    x,
    y,
    z,
    field=ux,
    field_name="Displacement ux (m)",
    case_dir=CASE_DIR,
    y_target=y_target,
    z_target=z_target,
    color="tab:green",
)
plot_field_along_x(
    x,
    y,
    z,
    field=sigma_vm,
    field_name=r"Von Mises stress (Pa)",
    case_dir=CASE_DIR,
    y_target=y_target,
    z_target=z_target,
    color="tab:orange",
)

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
sigma_xx_num = np.mean(sigma_xx)
sigma_vm_num = np.mean(sigma_vm)
ux_num = np.max(ux)

errors = {
    "sigma_xx": {
        "numerical": float(sigma_xx_num),
        "reference": SIGMA_REF,
        "abs_error": float(abs(sigma_xx_num - SIGMA_REF)),
        "rel_error": float(abs(sigma_xx_num - SIGMA_REF) / SIGMA_REF),
    },
    "sigma_von_mises": {
        "numerical": float(sigma_vm_num),
        "reference": VON_MISES_REF,
        "abs_error": float(abs(sigma_vm_num - VON_MISES_REF)),
        "rel_error": float(abs(sigma_vm_num - VON_MISES_REF) / VON_MISES_REF),
    },
    "ux_displacement": {
        "numerical": float(ux_num),
        "reference": UX_REF,
        "abs_error": float(abs(ux_num - UX_REF)),
        "rel_error": float(abs(ux_num - UX_REF) / UX_REF),
    },
}

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] non-regression completed.\n")
