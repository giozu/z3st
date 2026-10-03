#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: spherical_shell

Steady thermo-elasticity of a spherical shell (Ri = 2.0 m, Ro = 2.5 m) heated
by gamma radiation, q(r) = q0 (Ri/r) exp(-mu (r - Ri)), with fixed temperatures
on both surfaces, a free inner surface and a clamped outer surface. One octant
is meshed with structured hexahedra, graded radially towards the inner surface.

Compared against the exact solution: T(r) in closed form, and the
thermo-elastic displacement and stresses of a sphere for that T(r).

non-regression script
---------------------
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

from z3st.utils.non_regression import case_paths, error_metric, finish, load_case, metric
from z3st.utils.utils_extract_vtu import extract_principal_stresses, extract_temperature, list_fields
from z3st.utils.utils_plot import plot_field_along_r_xyz

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR, VTU_FILE, OUT_JSON = case_paths(__file__)

# Parameters
geom, inp, mat = load_case(CASE_DIR)
Ri = float(geom["Ri"])  # inner radius (m)
Ro = float(geom["Ro"])  # outer radius (m)
t = Ro - Ri  # thickness (m)


Pi = 0.000  # internal pressure (Pa)
Po = 0.000  # external pressure (Pa)

nu = float(mat["nu"])  # Poisson's ratio
E = float(mat["E"])  # Young modulus (Pa)

alpha = float(mat["alpha"])  # (1/K)
k = float(mat["k"])  # (W/mK
# Exact solution and every parameter it needs: exact.py in this directory.
sys.path.insert(0, CASE_DIR)
import exact  # noqa: E402

P = exact.parameters(CASE_DIR)
Ti, T0 = P["Ti"], P["To"]

# Default mesh (8 x 8 x 20 hexahedra): largest error 8.8e-3 (Linf of T over the
# temperature rise). Refinement series 4/8/16 cells per patch edge and 10/20/40
# radial layers: T, sigma_rr, sigma_tt and u_r all converge at second order.
TOLERANCE = 2e-2

# ANSI colors
GREEN = "\033[92m"
RED = "\033[91m"
BOLD = "\033[1m"
END = "\033[0m"


# --.. ..- .-.. .-.. --- analytical solution --.. ..- .-.. .-.. ---
def analytic_T(r):
    return exact.temperature(r, P)


def exact_mechanics(r_eval):
    return exact.mechanics(r_eval, P)


def sigma_th(x, T_num):
    """Approximate analytical thermal stress from the temperature profile."""
    # Constraints in 1, 2 or 3 directions:
    # c = 0, 1, 2, respectively
    c = 3.0

    T_mean = 3.0 / (Ro**3 - Ri**3) * np.trapezoid(T_num * x**2, x)
    return alpha * E / (1.0 - c * nu) * (T_mean - T_num)


# --.. ..- .-.. .-.. --- checks --.. ..- .-.. .-.. ---
if not os.path.exists(VTU_FILE):
    raise FileNotFoundError(f"[ERROR] VTU file not found: {VTU_FILE}")
print(f"[INFO] Using VTU file: {VTU_FILE}")

# --.. ..- .-.. .-.. --- list fields --.. ..- .-.. .-.. ---
list_fields(VTU_FILE)

# --. mesh bounds --..
grid = pv.read(VTU_FILE)
xmin, xmax, ymin, ymax, zmin, zmax = grid.bounds
print(f"\n[INFO] Grid bounds:")
print(f"\tx ∈ [{xmin:.4e}, {xmax:.4e}]")
print(f"\ty ∈ [{ymin:.4e}, {ymax:.4e}]")
print(f"\tz ∈ [{zmin:.4e}, {zmax:.4e}]")

# --.. ..- .-.. .-.. --- extract field --.. ..- .-.. .-.. ---
x, y, z, T = extract_temperature(VTU_FILE)

r_ref = np.linspace(Ri, Ro, 400)
T_ref = analytic_T(r_ref)
slenderness = 0.5 * (Ri + Ro) / t

rr, TT = plot_field_along_r_xyz(
    x,
    y,
    z,
    T,
    f"Temperature (K), R/t = {slenderness:.2f}",
    CASE_DIR,
    color="#0072B2",
    average="round",  # "round" or False
    decimals=2,
    r_ref=r_ref,
    f_ref=T_ref,
    label_ref="Analytical solution (sphere)",
)

r, _, _, sigma1, sigma2, sigma3 = extract_principal_stresses(
    VTU_FILE,
    average="round",  # "round" or None
    decimals=2,
    return_coords=True,
)

# --.. ..- .-.. .-.. --- analytical comparison --.. ..- .-.. .-.. ---
plt.figure()

plt.plot(r, sigma1, "-", color="#D55E00", lw=2, label=r"$\sigma_1$")
plt.plot(r, sigma2, "-", color="#009E73", lw=2, label=r"$\sigma_2$")
plt.plot(r, sigma3, "-", color="#0072B2", lw=2, label=r"$\sigma_3$")

# thermal_stress = sigma_th(rr, TT)
# plt.plot(rr, thermal_stress, color="#CC79A7", lw=2, label=r"$\sigma_{th}$")

plt.xlabel("r (m)")
plt.ylabel("Stress (Pa)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(CASE_DIR, "output", "stress_comparison.png"), dpi=200)
plt.close()

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
Tmax_num = float(np.max(T))
Tmax_ref = float(np.max(T_ref))
# pointwise error over every node of the shell
r_nodes = np.sqrt(x**2 + y**2 + z**2)
dT = T - analytic_T(r_nodes)
L2_T = float(np.sqrt(np.mean(dT**2)))
RelL2_T = L2_T / float(np.sqrt(np.mean(analytic_T(r_nodes) ** 2)))
Linf_T = float(np.max(np.abs(dT)))
rise = Tmax_ref - min(Ti, T0)
print(f"[INFO] T_max: numerical {Tmax_num:.3f} K, exact {Tmax_ref:.3f} K")
print(f"[INFO] max |T - T_exact| = {Linf_T:.3f} K, {Linf_T / rise:.2%} of the temperature rise")
RelErr_Tmax = abs(Tmax_num - Tmax_ref) / (Tmax_ref if Tmax_ref != 0 else 1.0)

# pointwise comparison of displacement (nodes) and stresses (cell centres)
cells = grid.cell_centers()
rc = np.linalg.norm(cells.points, axis=1)
nc = cells.points / rc[:, None]
S = np.asarray(grid.cell_data["Stress (cells)"]).reshape(-1, 3, 3)
srr_num = np.einsum("ci,cij,cj->c", nc, S, nc)
stt_num = 0.5 * (np.einsum("cii->c", S) - srr_num)
_, srr_ex, stt_ex = exact_mechanics(rc)
s_peak = float(np.max(np.abs(stt_ex)))
L2_srr = float(np.sqrt(np.mean((srr_num - srr_ex) ** 2)))
L2_stt = float(np.sqrt(np.mean((stt_num - stt_ex) ** 2)))

pn = grid.points
rn = np.linalg.norm(pn, axis=1)
ur_num = np.einsum("ni,ni->n", np.asarray(grid.point_data["Displacement"]), pn / rn[:, None])
ur_ex, _, _ = exact_mechanics(rn)
u_peak = float(np.max(np.abs(ur_ex)))
L2_ur = float(np.sqrt(np.mean((ur_num - ur_ex) ** 2)))
print(f"[INFO] peak |sigma_tt| = {s_peak / 1e6:.3f} MPa, peak |u_r| = {u_peak * 1e6:.3f} um")
print(f"[INFO] L2 error / peak: sigma_rr {L2_srr / s_peak:.2e}, sigma_tt {L2_stt / s_peak:.2e}, u_r {L2_ur / u_peak:.2e}")

errors = {
    "T_max": metric(Tmax_num, Tmax_ref, rel=RelErr_Tmax),
    "L2_error_T": error_metric(L2_T, rel=RelL2_T),
    "Linf_error_T": error_metric(Linf_T, rel=Linf_T / rise),
    "L2_error_sigma_rr": error_metric(L2_srr, rel=L2_srr / s_peak),
    "L2_error_sigma_tt": error_metric(L2_stt, rel=L2_stt / s_peak),
    "L2_error_u_r": error_metric(L2_ur, rel=L2_ur / u_peak),
}


finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)
