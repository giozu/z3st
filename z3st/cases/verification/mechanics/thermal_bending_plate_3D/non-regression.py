#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: verification/mechanics/thermal_bending_plate_3D

non-regression script
---------------------
Free flat plate (all faces traction-free) with a linear temperature
T(x) = Ti + g x across the thickness, g = (To - Ti) / Lx.

The free-plate stress (notebook "thermal stress in flat plate.ipynb")

    sigma_yy = sigma_zz = alpha E / (1 - nu) * [ T_mean
                          + 12 (x - L/2) / L^3 * int_0^L T (x - L/2) dx - T(x) ]

is identically zero for a linear T: the plate takes the gradient up by
bending, with membrane strain eps_0 = alpha (T_mean - T_ref) and curvature
kappa = alpha g. The exact displacement, up to an x-translation, is

    u_x = alpha [ (Ti - T_ref) x + g x^2 / 2 - g (y^2 + z^2) / 2 ]
    u_y = alpha (T(x) - T_ref) y
    u_z = alpha (T(x) - T_ref) z

and it holds over the whole plate, edges included.
"""

import numpy as np

from z3st.utils.non_regression import case_paths, error_metric, finish, load_case, metric
from z3st.utils.utils_extract_vtu import extract_displacement, extract_stress, extract_temperature, list_fields

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR, VTU_FILE, OUT_JSON = case_paths(__file__)

geom, inp, mat = load_case(CASE_DIR)
Lx, Ly = float(geom["Lx"]), float(geom["Ly"])  # m   thickness, half-span

E, nu, alpha, T_ref = float(mat["E"]), float(mat["nu"]), float(mat["alpha"]), float(mat["T_ref"])  # Pa, -, 1/K, K
Ti, To = 600.0, 400.0  # K   face temperatures at x = 0 and x = Lx
g = (To - Ti) / Lx  # K/m temperature gradient

kappa_ref = alpha * g  # 1/m
eps0_ref = alpha * (0.5 * (Ti + To) - T_ref)  # -
sigma_scale = alpha * E / (1.0 - nu) * abs(To - Ti)  # Pa  stress if bending and expansion were both restrained

TOLERANCE = 2e-2  # -


# --.. ..- .-.. .-.. --- analytic functions  --.. ..- .-.. .-.. ---
def analytic_T(x):
    return Ti + g * x


def analytic_u(x, y, z):
    dT = analytic_T(x) - T_ref
    ux = alpha * ((Ti - T_ref) * x + 0.5 * g * x**2 - 0.5 * g * (y**2 + z**2))
    return np.column_stack([ux, alpha * dT * y, alpha * dT * z])


# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
list_fields(VTU_FILE)

x_T, _, _, T = extract_temperature(VTU_FILE)
T_ref_x = analytic_T(x_T)

x, y, z, u = extract_displacement(VTU_FILE)
u_ref = analytic_u(x, y, z)
u_ref[:, 0] += np.mean(u[:, 0] - u_ref[:, 0])  # x-translation is not fixed by the BCs

# Curvature and membrane strain on the mid-plane line x = Lx/2, z = 0
line = (np.abs(x - Lx / 2) < 1e-9) & (np.abs(z) < 1e-9)
kappa_num = -2.0 * np.polyfit(y[line], u[line, 0], 2)[0]
tip = line & (np.abs(y - Ly) < 1e-9)
eps0_num = float(u[tip, 1][0] / Ly)

_, _, _, s = extract_stress(VTU_FILE, component="all", return_coords=True, prefer="cells")
sigma_rms = {c: float(np.sqrt(np.mean(s[c] ** 2))) for c in ("xx", "yy", "zz")}

print(f"[INFO] kappa: numerical {kappa_num:.6e} 1/m, analytic {kappa_ref:.6e} 1/m")
print(f"[INFO] eps_0: numerical {eps0_num:.6e}, analytic {eps0_ref:.6e}")

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
L2_T = float(np.sqrt(np.mean((T - T_ref_x) ** 2)))
L2_u = float(np.sqrt(np.mean(np.sum((u - u_ref) ** 2, axis=1))) / np.max(np.linalg.norm(u_ref, axis=1)))

errors = {
    "L2_error_T": error_metric(L2_T, rel=L2_T / np.mean(np.abs(T_ref_x))),
    "kappa": metric(kappa_num, kappa_ref, rel=abs(kappa_num - kappa_ref) / abs(kappa_ref)),
    "eps_0": metric(eps0_num, eps0_ref, rel=abs(eps0_num - eps0_ref) / abs(eps0_ref)),
    "L2_error_u": error_metric(L2_u),
    **{f"rms_sigma_{c}": error_metric(v, rel=v / sigma_scale) for c, v in sigma_rms.items()},
}

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)
