#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: U_coaxial_contact_2D (PCMI penalty contact)

non-regression
--------------
Compares the Z3ST penalty contact pressure against the analytical Lame
interference-fit pressure over the power ramp.

The interference delta = u_pellet_free(b) - u_clad_free(b) - g0 is estimated
from the thermal free expansion (plane-strain, area-mean temperature); the
Lame shrink-fit relation then gives the contact pressure. This is an
*approximate* reference for the thermal case (the interference is a small
difference of large expansions); a tight verification uses an isothermal case
with a prescribed geometric interference.
"""

import os
import glob
import yaml
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv

from z3st.utils.utils_load import generate_power_history

CASE = os.path.dirname(__file__)
OUT = os.path.join(CASE, "output")
files = sorted(glob.glob(os.path.join(OUT, "fields_*.vtu")))

# --- geometry, contact and material parameters ---
geo = yaml.safe_load(open(os.path.join(CASE, "geometry.yaml")))
b = float(geo["outer_radius_1"])     # interface radius (pellet outer ~ clad inner)
bci = float(geo["inner_radius_2"])
c = float(geo["outer_radius_2"])

inp = yaml.safe_load(open(os.path.join(CASE, "input.yaml")))
cc = inp["models"]["contact"]
g0 = float(cc["initial_gap"])
k_pen = float(cc["penalty_stiffness"])
fuel = yaml.safe_load(open(os.path.join(CASE, inp["materials"]["cyl_1"])))
clad = yaml.safe_load(open(os.path.join(CASE, inp["materials"]["cyl_2"])))
Ef, nuf, af, Trf = (float(fuel[k]) for k in ("E", "nu", "alpha", "T_ref"))
Ec, nuc, ac, Trc = (float(clad[k]) for k in ("E", "nu", "alpha", "T_ref"))

times, _, _ = generate_power_history(inp["time"], inp["lhr"],
                                     n_steps=inp["n_steps"] - 1, filename=None)

# plane-strain effective elastic properties
ps = lambda E, nu: (E / (1 - nu**2), nu / (1 - nu))
Efs, nufs = ps(Ef, nuf)
Ecs, nucs = ps(Ec, nuc)
# Lame interference-fit compliance: delta = p * b * comp  (solid shaft in tube)
comp = (1.0 / Ecs) * ((c**2 + bci**2) / (c**2 - bci**2) + nucs) + (1.0 / Efs) * (1.0 - nufs)

# --- per-step contact pressure: Z3ST penalty vs analytical Lame ---
surf = lambda v, r, r0: v[np.abs(r - r0) < 2e-5].mean()      # surface mean
amean = lambda v, r, lo, hi: (lambda m: np.sum(v[m] * r[m]) / np.sum(r[m]))((r >= lo) & (r <= hi))  # area mean

p_z3st, p_lame = [], []
for f in files:
    m = pv.read(f)
    r = m.points[:, 0]
    ur = m.point_data["Displacement"][:, 0]
    T = m.point_data["Temperature"]

    gap = g0 + surf(ur, r, bci) - surf(ur, r, b)            # current gap
    p_z3st.append(k_pen * max(0.0, -gap) / 1e6)             # MPa

    # free thermal interference (plane-strain, area-mean temperature)
    u_f = af * (1 + nuf) * b * (amean(T, r, 0.0, b) - Trf)
    u_c = ac * (1 + nuc) * bci * (amean(T, r, bci, c) - Trc)
    delta = u_f - u_c - g0
    p_lame.append((delta / (b * comp) / 1e6) if delta > 0 else 0.0)

p_z3st, p_lame, t = np.array(p_z3st), np.array(p_lame), np.array(times)

# --- single plot: contact pressure over time ---
plt.figure(figsize=(7, 5))
plt.plot(t, p_lame, "k--", lw=1.5, label="Analytical Lame interference (approx.)")
plt.plot(t, p_z3st, "r-o", lw=2, label="Z3ST penalty contact")
plt.xlabel("time (s)")
plt.ylabel("contact pressure (MPa)")
plt.title("Contact pressure: Z3ST vs analytical Lame")
plt.grid(True, ls=":", alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT, "contact_pressure_verification.png"), dpi=150)
print("[INFO] contact_pressure_verification.png saved")

mask = p_lame > 1.0
if mask.any():
    rel = np.abs(p_z3st[mask] - p_lame[mask]) / p_lame[mask]
    print(f"[INFO] closed-gap steps: {mask.sum()}, mean rel. diff = {rel.mean() * 100:.1f}%")
print("[INFO] non-regression completed.\n")
