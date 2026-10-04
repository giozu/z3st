#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: teaching/03_plate_with_hole_2D  --  Kirsch: a circular hole in a
wide plate under remote tension.

Remote tension sigma along x, hole of radius a at the origin, a quarter of the
plate meshed with x = 0 and y = 0 as symmetry planes. Kirsch's solution for an
infinite plate (Timoshenko and Goodier, Theory of Elasticity, sec. 35), with
phi measured from the load direction:

    on the hole edge     sigma_tt(a, phi) = sigma (1 - 2 cos 2 phi)
                         -> 3 sigma at phi = 90 deg, -sigma at phi = 0
    across the section   sigma_xx(0, r) = sigma (1 + a^2/(2 r^2) + 3 a^4/(2 r^4))
                         -> 1.22 sigma at r = 2a, 1.02 sigma at r = 5a

The regime is plane strain (regime: 2d, eps_zz = 0). The in-plane stresses of
this problem do not depend on the elastic constants, so they are Kirsch's in
plane strain and in plane stress alike; what plane strain adds is
sigma_zz = nu (sigma_xx + sigma_yy), and through it the von Mises stress.

The plate is finite (half-side 20 a) and the mesh is too: the tolerance is
3 %, and the errors actually reached are printed.
"""

import numpy as np

from z3st.utils.non_regression import case_paths, finish, load_case, load_yaml, metric

CASE_DIR, VTU_FILE, OUT_JSON = case_paths(__file__)
geom, inp, mat = load_case(CASE_DIR)
bcs = load_yaml(CASE_DIR, "boundary_conditions.yaml")["mechanical"]["plate"]

a = float(geom["a"])
nu = float(mat["nu"])
P = next(float(b["traction"]) for b in bcs if b["type"] == "Neumann" and b["region"] == "xmax")
TOLERANCE = 3e-2


def kirsch_section(r):
    """sigma_xx along x = 0, the section through the hole normal to the load."""
    return P * (1 + 0.5 * (a / r) ** 2 + 1.5 * (a / r) ** 4)


def kirsch_edge(phi):
    """Hoop stress on the hole edge, phi from the load direction."""
    return P * (1 - 2 * np.cos(2 * phi))


# --.. ..- .-.. .-.. --- FE fields at the nodes --.. ..- .-.. .-.. ---
import pyvista as pv  # noqa: E402

grid = pv.read(VTU_FILE)
x, y = grid.points[:, 0], grid.points[:, 1]
r = np.hypot(x, y)
S = np.asarray(grid.point_data["Stress (points)"]).reshape(-1, 3, 3)


def at(px, py):
    return int(np.argmin(np.hypot(x - px, y - py)))


i_top, i_side = at(0.0, a), at(a, 0.0)
Kt = S[i_top, 0, 0] / P
sigma_side = S[i_side, 1, 1] / P

# the section x = 0, from the hole edge to 5 a
sec = (np.abs(x) < 1e-9) & (r <= 5 * a + 1e-9)
order = np.argsort(r[sec])
r_sec, sxx_sec = r[sec][order], S[sec, 0, 0][order]
sec_err = np.sqrt(np.mean((sxx_sec - kirsch_section(r_sec)) ** 2)) / P

# the hole edge
edge = np.abs(r - a) < 1e-9
phi = np.arctan2(y[edge], x[edge])
et = np.c_[-np.sin(phi), np.cos(phi), np.zeros_like(phi)]
stt = np.einsum("ni,nij,nj->n", et, S[edge], et)
edge_err = np.sqrt(np.mean((stt - kirsch_edge(phi)) ** 2)) / P

print(f"\nhole a = {a*1e3:.0f} mm, remote tension sigma = {P/1e6:.0f} MPa")
print(f"\non the hole edge, 90 deg from the load : sigma_xx = {Kt:.3f} sigma   (Kirsch: 3)")
print(f"on the hole edge, on the load axis     : sigma_yy = {sigma_side:+.3f} sigma  (Kirsch: -1)")
print("\nacross the section x = 0:")
for rr in (1, 2, 3, 5):
    j = int(np.argmin(np.abs(r_sec - rr * a)))
    print(f"   r = {r_sec[j]/a:4.2f} a   sigma_xx = {sxx_sec[j]/P:5.3f} sigma"
          f"   (Kirsch {kirsch_section(r_sec[j])/P:5.3f})")
print(f"\nRMS error / sigma: section {sec_err:.2e}, hole edge {edge_err:.2e}")
print(f"plane strain: sigma_zz at the peak = {S[i_top, 2, 2]/P:.3f} sigma"
      f" = nu (sigma_xx + sigma_yy) = {nu*(S[i_top, 0, 0] + S[i_top, 1, 1])/P:.3f} sigma")

errors = {
    "Kt": metric(Kt, 3.0),
    "sigma_on_load_axis": metric(sigma_side, -1.0),
    "section_profile_rms": {"numerical": sec_err, "reference": 0.0,
                            "abs_error": sec_err, "rel_error": sec_err},
    "hole_edge_rms": {"numerical": edge_err, "reference": 0.0,
                      "abs_error": edge_err, "rel_error": edge_err},
}

finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)
