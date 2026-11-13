#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: hydrostatic_cylinder

Open thin-walled cylinder filled with a heavy liquid (Neumann BC with
`hydrostatic:`). The wall pressure is p(z) = rho*g*(z_free - z). With no
external pressure, radial equilibrium integrated through the wall gives

    int sigma_tt dr = p(z) * Ri + d/dz int r sigma_rz dr,

whatever the thickness. With p linear in z the membrane solution needs no
bending, so the shear term vanishes and N_theta = p * Ri, except in a
bending layer at the base, where Clamp_y stops the wall from rotating.
The L2 norm over the whole height keeps that layer small.
"""

import os

import matplotlib.pyplot as plt
import numpy as np

from z3st.utils.non_regression import case_paths, error_metric, finish, load_case
from z3st.utils.utils_extract_vtu import extract_field

CASE_DIR, VTU_FILE, OUT_JSON = case_paths(__file__)
geom, inp, mat = load_case(CASE_DIR)
Ri, Ro = float(geom["Ri"]), float(geom["Ro"])
t = Ro - Ri

import yaml
with open(os.path.join(CASE_DIR, "boundary_conditions.yaml")) as f:
    hydro = next(bc["hydrostatic"] for bc in yaml.safe_load(f)["mechanical"]["steel"] if "hydrostatic" in bc)
rho_f, z_free, g = float(hydro["rho"]), float(hydro["z_free"]), 9.81

TOLERANCE = 1.0e-2

# Hoop force per unit length at each cell row (structured mesh: uniform dr).
x, z, _, S = extract_field(VTU_FILE, field_name="Stress (cells)")
z_rows = np.unique(np.round(z, 9))
N_fe = np.array([S[np.isclose(z, zr), 4].mean() * t for zr in z_rows])
N_ana = rho_f * g * (z_free - z_rows) * Ri
# Clamp_y holds u_z on the whole bottom edge, so it stops the wall from rotating
# there: a local bending layer, about 3*sqrt(Ri*t) high, carries sigma_zz.
z_edge = 3.0 * np.sqrt(Ri * t)
sigma_zz_max = np.abs(S[z > z_edge, 8]).max()

err_N = np.sqrt(np.mean((N_fe - N_ana) ** 2)) / np.sqrt(np.mean(N_ana**2))
err_zz = sigma_zz_max / (rho_f * g * z_free * Ri / t)  # axial stress vanishes (open top) above the edge layer
print(f"[INFO] relative L2 error on N_theta(z): {err_N:.3e}")
print(f"[INFO] max |sigma_zz| / max hoop:        {err_zz:.3e}")

plt.figure(figsize=(7, 5))
plt.plot(N_ana / t * 1e-6, z_rows, "-", color="#0072B2", label=r"$\rho g (z_f - z) R_i / t$")
plt.plot(N_fe / t * 1e-6, z_rows, "o", color="#D55E00", ms=3, label="Z3ST, through-wall average")
plt.xlabel(r"Average hoop stress $\bar\sigma_{\theta\theta}$ (MPa)")
plt.ylabel("z (m)")
plt.grid(True, linestyle="--", alpha=0.7)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(CASE_DIR, "output", "hoop_vs_height.png"), dpi=200)

errors = {
    "L2_error_N_theta": error_metric(err_N),
    "max_sigma_zz_over_hoop": error_metric(err_zz),
}
finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)
