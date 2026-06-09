#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: V_fuel_swelling_verification

Verifies that the eigenstrain bus consumes the burnup state bus — i.e. a fuel
material whose swelling is a function of the accumulated burnup field. A fissile
block accumulates a uniform burnup (flat power profile) over a fixed time; its
swelling eigenstrain ΔV/V = swelling_rate * bu then drives a free, stress-free
expansion. With alpha absent (no thermal strain) the closed form is

    bu      = q''' * t / (rho * HM * 8.64e10)        [MWd/kgU]
    ΔV/V    = swelling_rate * bu
    u_x(Lx) = (ΔV/V)/3 * Lx     and     sigma ≈ 0

q''' = lhr / area, area = Lx*Ly, 8.64e10 = 86400 s/day * 1e6 W/MW. This is the
chain state -> eigenstrain -> stress: it confirms swelling reads the burnup
field and produces the right strain, the analogue of V_swelling_verification but
for a *state-dependent* eigenstrain.
"""

import os
import glob
import yaml
import numpy as np

from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_verification import *

CASE_DIR = os.path.dirname(__file__)
OUT = os.path.join(CASE_DIR, "output")
OUT_JSON = os.path.join(OUT, "non-regression.json")

_single = os.path.join(OUT, "fields.vtu")
_steps = sorted(glob.glob(os.path.join(OUT, "fields_*.vtu")))
VTU_FILE = _steps[-1] if _steps else _single

# --. geometry + material + history --..
with open(os.path.join(CASE_DIR, "geometry.yaml")) as f:
    geom = yaml.safe_load(f)
Lx, Ly = float(geom["Lx"]), float(geom["Ly"])
area = Lx * Ly

with open(os.path.join(CASE_DIR, "input.yaml")) as f:
    inp = yaml.safe_load(f)
lhr = float(inp["lhr"][-1])
t_total = float(max(inp["time"]))

with open(os.path.join(CASE_DIR, next(iter(inp["materials"].values())))) as f:
    mat = yaml.safe_load(f)
rho = float(mat["rho"])
hm = float(mat.get("heavy_metal_fraction", 0.8815))
rate = float(mat["swelling_rate"])
E, nu = float(mat["E"]), float(mat["nu"])

# --. analytical references --..
SECONDS_PER_MWD = 8.64e10
q_vol = lhr / area
bu_ref = q_vol * t_total / (rho * hm * SECONDS_PER_MWD)   # uniform burnup (MWd/kgU)
dVV = rate * bu_ref                                       # volumetric swelling
eps_lin = dVV / 3.0                                       # linear eigenstrain / dir
UX_REF = eps_lin * Lx                                     # free expansion of xmax
K = E / (3.0 * (1.0 - 2.0 * nu))                          # bulk modulus
TOLERANCE = 1e-2

# --. numerical results --..
x_u, y_u, z_u, u = extract_field(VTU_FILE, field_name="Displacement")
ux_num = float(np.max(u[:, 0]))

_, _, _, vm = extract_field(VTU_FILE, field_name="VonMises (cells)")
vm_num = float(np.mean(np.abs(vm)))                       # ≈ 0 for free swelling

print(f"[INFO] burnup (uniform) = {bu_ref:.4f} MWd/kgU  ->  ΔV/V = {dVV:.4e}")
print(f"[INFO] u_x(max): numerical = {ux_num:.6e} m, analytical = {UX_REF:.6e} m")
print(f"[INFO] mean |von Mises|: {vm_num:.3e} Pa  (free swelling ⇒ ≈ 0)")
print(f"[INFO] fully constrained would give σ = -K·ΔV/V = {-K*dVV:.3e} Pa")

errors = {
    "u_x_free_swelling": {
        "numerical": ux_num,
        "reference": UX_REF,
        "abs_error": float(abs(ux_num - UX_REF)),
        "rel_error": float(abs(ux_num - UX_REF) / UX_REF),
    },
    "stress_free_residual": {
        "numerical": vm_num,
        "reference": 0.0,
        "abs_error": vm_num,
        "rel_error": float(vm_num / (K * dVV)),
    },
}

pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] non-regression completed.\n")
