#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: teaching/02_tensile_bar_3D  --  a round tensile bar, solved in 3D.

The specimen of a uniaxial tensile test: d = 10 mm, L = 50 mm, pulled by
F = 15 kN, steel (E = 200 GPa, nu = 0.3). The whole bar is meshed, held on its interior
planes x = 0 and y = 0 (which the exact solution leaves in place) and axially at
the bottom end, so nothing restrains the lateral contraction and the bar is in
uniaxial stress:

    sigma_zz = F/A,  every other component 0
    eps_zz = sigma/E,  eps_rr = -nu sigma/E,  eps_v = (1 - 2 nu) sigma/E
    on a plane whose normal is at theta from the axis:
        sigma_n = sigma cos^2 theta,   tau = (sigma/2) sin 2 theta
    Tresca = von Mises = sigma,  so yield starts at F_y = S_y A
    in a direction at theta from the axis (Mohr's circle of strain):
        eps_n = eps_zz cos^2 theta + eps_rr sin^2 theta
        gamma/2 = (eps_zz - eps_rr) sin theta cos theta,  largest at 45 deg,
        where gamma = tau_max / G

The displacement field is linear in (x, y, z), so linear tetrahedra reproduce
it to machine precision on any mesh. The numbers printed are the ones a hand
calculation gives; here they are read off the FE field instead. The figures
are drawn by plots.py.
"""

import numpy as np

from z3st.utils.non_regression import case_paths, finish, load_case, load_yaml, metric
from z3st.utils.utils_extract_vtu import extract_displacement, extract_field

CASE_DIR, VTU_FILE, OUT_JSON = case_paths(__file__)
geom, inp, mat = load_case(CASE_DIR)
bcs = load_yaml(CASE_DIR, "boundary_conditions.yaml")["mechanical"]["steel"]

R, L = float(geom["Ro"]), float(geom["Lz"])
E, nu = float(mat["E"]), float(mat["nu"])
P = next(float(b["traction"]) for b in bcs if b["type"] == "Neumann" and b["region"] == "top")
S_y = 235e6  # Pa, S235: the 235 is the grade

A = np.pi * R**2
F = P * A
TOLERANCE = 1e-6

# --.. ..- .-.. .-.. --- analytical reference --.. ..- .-.. .-.. ---
sigma_ref = P
eps_z_ref = P / E
eps_r_ref = -nu * P / E

# --.. ..- .-.. .-.. --- FE fields --.. ..- .-.. .-.. ---
_, _, _, S = extract_field(VTU_FILE, field_name="Stress (cells)")
_, _, _, vm = extract_field(VTU_FILE, field_name="VonMises (cells)")
S = S.reshape(-1, 3, 3)
sig = S.mean(axis=0)                    # the state is uniform: one tensor describes it
_, _, _, eps_cells = extract_field(VTU_FILE, field_name="Strain (cells)")
eps = eps_cells.reshape(-1, 3, 3).mean(axis=0)   # the tensor strain, not engineering gamma
G = E / (2 * (1 + nu))

x, y, z, u = extract_displacement(VTU_FILE)
r = np.hypot(x, y)
top = np.abs(z - L) < 1e-9
rim = np.abs(r - R) < 1e-9
dL = float(u[top, 2].mean())
u_r = float(((x[rim] * u[rim, 0] + y[rim] * u[rim, 1]) / r[rim]).mean())
eps_z, eps_r = dL / L, u_r / R
eps_v = eps_z + 2 * eps_r

w = np.linalg.eigvalsh(sig)
tresca = w[-1] - w[0]
von_mises = float(vm.mean())

# --.. ..- .-.. .-.. --- the hand calculation, from the FE field --.. ..- .-.. .-.. ---
print(f"\nbar: d = {2*R*1e3:.0f} mm, L = {L*1e3:.0f} mm, A = {A*1e6:.1f} mm^2, F = {F/1e3:.1f} kN")
print(f"\n1. stress        sigma_zz = {sig[2, 2]/1e6:.1f} MPa   (F/A = {P/1e6:.1f} MPa)")
print(f"   largest other component  {np.abs(sig - np.diag([0, 0, sig[2, 2]])).max()/1e6:.1e} MPa")

print("\n2. a plane whose normal is at theta from the axis: t = sigma^T n")
print(f"   {'theta':>6} {'sigma_n (MPa)':>14} {'tau (MPa)':>10}")
thetas = np.radians(np.linspace(0, 90, 181))
sn_fe, tau_fe = [], []
for th in thetas:
    n = np.array([np.sin(th), 0.0, np.cos(th)])
    t = sig.T @ n
    sn = t @ n
    sn_fe.append(sn)
    tau_fe.append(np.sqrt(max(t @ t - sn**2, 0.0)))
sn_fe, tau_fe = np.array(sn_fe), np.array(tau_fe)
for deg in (0, 30, 45, 60, 90):
    i = 2 * deg
    print(f"   {deg:>5}d {sn_fe[i]/1e6:14.1f} {tau_fe[i]/1e6:10.1f}")
print(f"   largest shear {tau_fe.max()/1e6:.1f} MPa at {np.degrees(thetas[tau_fe.argmax()]):.0f} deg: sigma/2")

print("\n3. strains, from the strain field and from the displacements of the ends:")
print(f"                 eps_zz = {eps[2, 2]:.3e}  ({eps_z:.3e})  ->  dL = {dL*1e3:.4f} mm")
print(f"                 eps_xx = eps_yy = {eps[0, 0]:.3e}  ({eps_r:.3e})  ->  dd = {2*u_r*1e3:.4f} mm")
print(f"                 eps_v  = {np.trace(eps):.3e}  = (1 - 2 nu) sigma/E: the volume grows")
gamma_max = eps[2, 2] - eps[0, 0]   # engineering shear at 45 deg, twice the tensor one
print(f"                 at 45 deg: gamma = eps_zz - eps_xx = {gamma_max:.3e}"
      f" = tau_max/G = {tau_fe.max()/G:.3e}")

print(f"\n4. yield         Tresca = {tresca/1e6:.1f} MPa,  von Mises = {von_mises/1e6:.1f} MPa")
print(f"                 both equal sigma in tension: yield at F_y = S_y A = {S_y*A/1e3:.1f} kN"
      f" (S_y = {S_y/1e6:.0f} MPa)")

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
errors = {
    "sigma_zz": metric(sig[2, 2], sigma_ref),
    "von_mises": metric(von_mises, sigma_ref),
    "tresca": metric(tresca, sigma_ref),
    "u_z_top": metric(dL, eps_z_ref * L),
    "u_r_rim": metric(u_r, eps_r_ref * R),
    "tau_max": metric(tau_fe.max(), sigma_ref / 2),
    "eps_zz": metric(eps[2, 2], eps_z_ref),
    "eps_xx": metric(eps[0, 0], eps_r_ref),
    "gamma_45": metric(gamma_max, sigma_ref / 2 / G),
}

finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)
