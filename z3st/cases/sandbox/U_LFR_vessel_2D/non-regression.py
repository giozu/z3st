#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST post-processing --.. ..- .-.. .-.. ---
"""
U_LFR_vessel_2D: LFR vessel (cylindrical shell + hemispherical bottom head)
hung from its top edge, loaded by the lead head and its own weight.

Hand values, all from equilibrium, on through-wall AVERAGES (so the
local bending at the junction and at the top edge drops out):
  - shell hoop:   int sigma_tt dr = p(z) Ri                    (ring)
  - shell axial:  sigma_zz_avg = (W_Pb + W_steel,below) / A     (everything
                  below a section hangs from it; the lead pressure on the
                  walls adds up to the whole lead weight)
  - head pole:    sigma_avg = p_pole Ri^2 / (Ro^2 - Ri^2)       (cap at the pole)
Bending is reported, not checked, and so is the shell hoop stress: the
shell is only ~7.5 decay lengths long, so the junction disturbance never
dies out. Its source is the head equator, where the hemisphere carries all
the hanging lead in meridional tension and Laplace, N_t = p R - N_phi,
turns its hoop stress compressive (about -10 MPa) against +35 MPa in the
cylinder.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import yaml

from z3st.utils.non_regression import case_paths, finish, load_case, metric, tracked
from z3st.utils.utils_extract_vtu import extract_field

CASE_DIR, VTU_FILE, OUT_JSON = case_paths(__file__)
geom, inp, mat = load_case(CASE_DIR)
Ri, Ro, H = float(geom["Ri"]), float(geom["Ro"]), float(geom["H"])
t, A = Ro - Ri, np.pi * (Ro**2 - Ri**2)
g = float(inp["mechanical"]["gravity"])
rho_s = float(mat["rho"])
with open(os.path.join(CASE_DIR, "boundary_conditions.yaml")) as f:
    hydro = next(bc["hydrostatic"] for bc in yaml.safe_load(f)["mechanical"]["steel"] if "hydrostatic" in bc)
rho_f, z_f = float(hydro["rho"]), float(hydro["z_free"])

p = lambda z: rho_f * g * np.maximum(z_f - z, 0.0)
V_Pb = 2 / 3 * np.pi * Ri**3 + np.pi * Ri**2 * (z_f - Ri)
W_Pb = rho_f * g * V_Pb
W_steel_below = lambda z: rho_s * g * (2 / 3 * np.pi * (Ro**3 - Ri**3) + A * (z - Ri))

MPa = 1e-6
TOLERANCE = 2.0e-2

x, z, _, S = extract_field(VTU_FILE, field_name="Stress (cells)")

# --- shell: through-wall averages row by row (structured mesh, uniform dr)
shell = z > Ri
z_rows = np.unique(np.round(z[shell], 9))
avg = lambda comp: np.array([S[shell & np.isclose(z, zr), comp].mean() for zr in z_rows])
s_tt, s_zz = avg(4), avg(8)
s_tt_ana = p(z_rows) * Ri / t
s_zz_ana = (W_Pb + W_steel_below(z_rows)) / A

# Bending: max |sigma_zz - average| over the wall at each row (surface minus membrane)
bend = np.array([np.abs(S[shell & np.isclose(z, zr), 8] - m).max() for zr, m in zip(z_rows, s_zz)])

# --- head pole: cells on the axis
pole = (x < 0.05) & (z < 0.0)
s_pole = 0.5 * (S[pole, 0].mean() + S[pole, 4].mean())
s_pole_ana = p(0.0) * Ri**2 / (Ro**2 - Ri**2)

# Mid-shell row: furthest from the junction and from the top edge
i_mid = np.argmin(np.abs(z_rows - 0.5 * (Ri + H)))

print("\n  Quantity                         hand (MPa)   Z3ST (MPa)")
print(f"  pole, membrane                   {s_pole_ana*MPa:10.2f}   {s_pole*MPa:10.2f}")
print(f"  shell hoop, z = {z_rows[i_mid]:.2f} m          {s_tt_ana[i_mid]*MPa:10.2f}   {s_tt[i_mid]*MPa:10.2f}")
print(f"  shell axial, z = {z_rows[i_mid]:.2f} m         {s_zz_ana[i_mid]*MPa:10.2f}   {s_zz[i_mid]*MPa:10.2f}")
print(f"  shell axial, top edge            {s_zz_ana[-1]*MPa:10.2f}   {s_zz[-1]*MPa:10.2f}")
s_eq_sph = p(Ri) * Ri / t - (W_Pb + W_steel_below(Ri)) / A
print(f"  head equator hoop (membrane)     {s_eq_sph*MPa:10.2f}   (cylinder side: {s_tt_ana[0]*MPa:.2f})")
print(f"  lead mass {rho_f*V_Pb/1e3:.0f} t, wall pressure at the pole {p(0.0)*MPa:.3f} MPa")
print(f"  max axial bending stress, shell  {bend.max()*MPa:10.2f}   at z = {z_rows[np.argmax(bend)]:.2f} m")
print(f"  max |sigma| anywhere             {np.abs(S).max()*MPa:10.2f}\n")

fig, ax = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
ax[0].plot(s_tt_ana * MPa, z_rows, "-", color="#0072B2", label="hand, p(z) Ri / t")
ax[0].plot(s_tt * MPa, z_rows, "o", color="#D55E00", ms=2, label="Z3ST, wall average")
ax[0].set_xlabel(r"Hoop stress $\bar\sigma_{\theta\theta}$ (MPa)")
ax[1].plot(s_zz_ana * MPa, z_rows, "-", color="#0072B2", label="hand, hanging weight / A")
ax[1].plot(s_zz * MPa, z_rows, "o", color="#D55E00", ms=2, label="Z3ST, wall average")
ax[1].plot(bend * MPa, z_rows, "--", color="#009E73", label="Z3ST, bending part")
ax[1].set_xlabel(r"Axial stress $\sigma_{zz}$ (MPa)")
ax[0].set_ylabel("z (m), shell")
for a in ax:
    a.grid(True, linestyle="--", alpha=0.7)
    a.legend(fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(CASE_DIR, "output", "shell_stresses.png"), dpi=200)

errors = {
    "pole_membrane": metric(s_pole, s_pole_ana),
    "shell_hoop_mid": tracked(s_tt[i_mid]),
    "shell_axial_mid": metric(s_zz[i_mid], s_zz_ana[i_mid]),
    "shell_axial_top": metric(s_zz[-1], s_zz_ana[-1]),
    "max_bending_shell": tracked(bend.max()),
}
finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)
