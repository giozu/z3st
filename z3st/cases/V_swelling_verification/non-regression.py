#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: V_swelling_verification

Verifies the eigenstrain bus. A block of a material
carrying a constant volumetric swelling ΔV/V is left free to expand from a
corner (symmetry pins on three faces). Free swelling is stress-free, so:

    σ = 0                         (no stress develops)
    u_x(x = Lx) = (ΔV/V)/3 · Lx   (uniform isotropic expansion)

since the isotropic eigenstrain is ε* = (ΔV/V)/3 · I (and tr(ε*) = ΔV/V). The
fully-constrained counterpart would give the hydrostatic σ = -K (ΔV/V), with
K the bulk modulus of the material.
"""

import os
import yaml
import numpy as np

from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_verification import *

CASE_DIR = os.path.dirname(__file__)
VTU_FILE = os.path.join(CASE_DIR, "output", "fields.vtu")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")

# --. geometry + material from the case YAML files --..
with open(os.path.join(CASE_DIR, "geometry.yaml")) as f:
    geom = yaml.safe_load(f)
Lx, Ly, Lz = float(geom["Lx"]), float(geom["Ly"]), float(geom["Lz"])
with open(os.path.join(CASE_DIR, "input.yaml")) as f:
    inp = yaml.safe_load(f)
with open(os.path.join(CASE_DIR, next(iter(inp["materials"].values())))) as f:
    mat = yaml.safe_load(f)
E, nu, s = float(mat["E"]), float(mat["nu"]), float(mat["swelling"])

# --. analytical references --..
eps_lin = s / 3.0                       # isotropic linear eigenstrain per direction
UX_REF = eps_lin * Lx                   # free expansion of the xmax face
K = E / (3.0 * (1.0 - 2.0 * nu))        # bulk modulus
TOLERANCE = 1e-2

# --. numerical results --..
x_u, y_u, z_u, u = extract_field(VTU_FILE, field_name="Displacement")
ux_num = float(np.max(u[:, 0]))

_, _, _, vm = extract_field(VTU_FILE, field_name="VonMises (cells)")
vm_num = float(np.mean(np.abs(vm)))     # should be close to zero for free swelling

print(f"[INFO] swelling ΔV/V = {s} → linear eigenstrain s/3 = {eps_lin:.4e}")
print(f"[INFO] u_x(max): numerical = {ux_num:.6e} m, analytical = {UX_REF:.6e} m")
print(f"[INFO] mean |von Mises|: {vm_num:.3e} Pa  (free swelling ⇒ ≈ 0)")
print(f"[INFO] bulk modulus K = {K:.3e} Pa; fully constrained would give σ = -K·s = {-K*s:.3e} Pa")

errors = {
    "u_x_free": {
        "numerical": ux_num,
        "reference": UX_REF,
        "abs_error": float(abs(ux_num - UX_REF)),
        "rel_error": float(abs(ux_num - UX_REF) / UX_REF),
    },
    "stress_free_residual": {
        # spurious stress normalised by the constrained-stress scale K·s
        "numerical": vm_num,
        "reference": 0.0,
        "abs_error": vm_num,
        "rel_error": float(vm_num / (K * s)),
    },
}

pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] non-regression completed.\n")
