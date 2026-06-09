#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: V_coaxial_contact_verification

VERIFICATION of the penalty contact pressure against the analytical Lame
interference-fit pressure. The pellet is heated UNIFORMLY, so its free thermal
expansion is exactly u(b) = alpha_f (T - T_ref) b (no gradient, no Poisson
ambiguity) and the cladding (held at T_ref) does not expand. The interference

    delta = alpha_f (T_pellet - T_ref) * b - g0

then gives the exact Lame shrink-fit pressure for a solid cylinder in a tube:

    p = delta / { b [ (1/E_c)((c^2+b^2)/(c^2-b^2) + nu_c) + (1/E_f)(1 - nu_f) ] }

(plane-stress form, consistent with the axially-free pellet). One plot:
Z3ST penalty contact pressure vs the Lame line, against pellet temperature.
"""

import os
import glob
import yaml
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv

from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_verification import *

CASE = os.path.dirname(__file__)
OUT = os.path.join(CASE, "output")
OUT_JSON = os.path.join(CASE, "output", "non-regression.json")

TOLERANCE = 5e-2

files = sorted(glob.glob(os.path.join(OUT, "fields_*.vtu")))

geo = yaml.safe_load(open(os.path.join(CASE, "geometry.yaml")))
b = float(geo["outer_radius_1"])     # pellet outer / interface radius
bci = float(geo["inner_radius_2"])
c = float(geo["outer_radius_2"])

inp = yaml.safe_load(open(os.path.join(CASE, "input.yaml")))
cc = inp["models"]["contact"]
g0 = float(cc["initial_gap"])
k_pen = float(cc["penalty_stiffness"])
fuel = yaml.safe_load(open(os.path.join(CASE, inp["materials"]["cyl_1"])))
Ef, nuf, af, Trf = (float(fuel[k]) for k in ("E", "nu", "alpha", "T_ref"))
clad = yaml.safe_load(open(os.path.join(CASE, inp["materials"]["cyl_2"])))
Ec, nuc = float(clad["E"]), float(clad["nu"])

# Lame interference-fit compliance (plane stress): delta = p * b * comp
comp = (1.0 / Ec) * ((c**2 + bci**2) / (c**2 - bci**2) + nuc) + (1.0 / Ef) * (1.0 - nuf)

surf = lambda v, r, r0: v[np.abs(r - r0) < 2e-5].mean()
amean = lambda v, r, lo, hi: (lambda m: np.sum(v[m] * r[m]) / np.sum(r[m]))((r >= lo) & (r <= hi))

T_pellet, p_z3st, p_lame = [], [], []
for f in files:
    m = pv.read(f)
    r = m.points[:, 0]
    ur = m.point_data["Displacement"][:, 0]
    T = m.point_data["Temperature"]

    Tp = amean(T, r, 0.0, b)                                  # uniform pellet temperature
    T_pellet.append(Tp)

    gap = g0 + surf(ur, r, bci) - surf(ur, r, b)
    p_z3st.append(k_pen * max(0.0, -gap) / 1e6)              # MPa

    delta = af * (Tp - Trf) * b - g0                          # exact interference
    p_lame.append((delta / (b * comp) / 1e6) if delta > 0 else 0.0)

T_pellet, p_z3st, p_lame = map(np.array, (T_pellet, p_z3st, p_lame))

plt.figure(figsize=(7, 5))
plt.plot(T_pellet, p_lame, "k--", lw=1.5, label="Analytical Lame interference (exact)")
plt.plot(T_pellet, p_z3st, "r-o", lw=2, label="Z3ST penalty contact")
plt.xlabel("pellet temperature (K)")
plt.ylabel("contact pressure (MPa)")
plt.title("Contact pressure verification: Z3ST vs Lame (uniform ΔT)")
plt.grid(True, ls=":", alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT, "contact_pressure_verification.png"), dpi=150)
print("[INFO] contact_pressure_verification.png saved")

mask = p_lame > 1.0
if mask.any():
    rel = np.abs(p_z3st[mask] - p_lame[mask]) / p_lame[mask]
    print(f"[INFO] closed-gap steps: {mask.sum()}, mean rel. error vs Lame = {rel.mean() * 100:.1f}%")
print("[INFO] non-regression completed.\n")

# --. numerical results --..
errors = {
    "contact_pressure": {
        "numerical": p_z3st[mask].max(),
        "reference": p_lame[mask].max(),
        "abs_error": float(rel.max()),
        "rel_error": float(rel.max()),
    },
}

pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE)
regression_check(errors, CASE)

print("\n[INFO] non-regression completed.\n")
