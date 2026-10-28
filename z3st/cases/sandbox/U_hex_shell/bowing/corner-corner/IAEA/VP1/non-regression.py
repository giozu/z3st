#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
# Author: Bianca Funaro
"""
Z3ST case: sandbox/U_hex_shell/bowing/corner-corner/IAEA/VP1

IAEA VP1 free thermal bow, constant 15-15Ti (alpha, E): one duct clamped at
the bottom, T linear corner to corner over the core and above it (Dirichlet
on every face, boundary_conditions.yaml). The hot corner is at -y, so the
duct bows toward +y. T is linear in y in every section, so the bow follows the
beam solution (verification/fuel_assembly/IAEA/vp1); the y·z term over the
core and the kinks at z = 1.5 and 2.5 m leave small stresses (~3 MPa):

    kappa(z) = alpha dT/dy r(z),   r(z) = clip(z - 1.5, 0, 1)
    u_y(z)   = int_0^z kappa(s) (z - s) ds
             = kappa0 (z - 1.5)^3 / 6                              core
             = kappa0 [1/6 + (z - 2.5)/2 + (z - 2.5)^2/2]           above

  - section-mean u_y at the core top, ACLP and duct top against the beam;
  - von Mises maximum and median: no closed form, tracked (gold only).
"""

import os

import matplotlib.pyplot as plt
import numpy as np

from z3st.utils.non_regression import finish, load_yaml, metric, tracked
from z3st.utils.utils_extract_xdmf import extract_field_xdmf

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(os.path.abspath(__file__))
inp = load_yaml(CASE_DIR, "input.yaml")
mat = load_yaml(CASE_DIR, inp["materials"]["steel"])
XDMF_FILE = os.path.join(CASE_DIR, "output", f"{inp['output']['filename']}.xdmf")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")

alpha = float(mat["alpha"])  # 1/K
dT_dy = 50.0 / 0.150         # K/m, 500 -> 550 °C over 0.150 m (boundary_conditions.yaml)
kappa0 = alpha * dT_dy       # 1/m, curvature above the core

STATIONS = {"core_top": 2.5, "ACLP": 3.0, "duct_top": 4.0}  # m

# z3st sits ~0.015 mm below the beam at every station, 1.4 % of the 1 mm at
# the core top (verification/fuel_assembly/IAEA/vp1: -2.1 % there).
TOLERANCE = 2.0e-2


# --.. ..- .-.. .-.. --- analytic bow --.. ..- .-.. .-.. ---
def beam_uy(z):
    z = np.asarray(z, dtype=float)
    c = np.clip(z - 1.5, 0.0, 1.0)
    a = np.clip(z - 2.5, 0.0, None)
    return kappa0 * (c**3 / 6.0 + a / 2.0 + a**2 / 2.0)


def von_mises(S):
    sxx, syy, szz = S[:, 0], S[:, 4], S[:, 8]
    sxy, syz, sxz = S[:, 1], S[:, 5], S[:, 2]
    return np.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2)
                   + 3.0 * (sxy**2 + syz**2 + sxz**2))


# --.. ..- .-.. .-.. --- read fields --.. ..- .-.. .-.. ---
_, _, zU, U = extract_field_xdmf(XDMF_FILE, field_name="Displacement", step_index=-1)
_, _, _, S = extract_field_xdmf(XDMF_FILE, field_name="Stress", step_index=-1)
vm = von_mises(S)

# Section mean of u_y at every node level
z_lev, inv = np.unique(np.round(zU, 9), return_inverse=True)
uy_mean = np.bincount(inv, weights=U[:, 1]) / np.bincount(inv)

errors = {}
for key, zs in STATIONS.items():
    num, ref = np.interp(zs, z_lev, uy_mean), float(beam_uy(zs))
    print(f"[INFO] u_y at {key:9s} (z = {zs:.1f} m): z3st = {num * 1e3:.4f} mm, beam = {ref * 1e3:.4f} mm")
    errors[f"uy_{key}"] = metric(num, ref)
print(f"[INFO] von Mises: max = {vm.max() * 1e-6:.3f} MPa, median = {np.median(vm) * 1e-6:.3f} MPa")
errors["vonMises_max"] = tracked(vm.max())
errors["vonMises_median"] = tracked(np.median(vm))

# --.. ..- .-.. .-.. --- plot --.. ..- .-.. .-.. ---
plt.figure(figsize=(6, 7))
plt.plot(uy_mean * 1e3, z_lev, "-", color="#0072B2", lw=2, label="z3st, section mean")
plt.plot(beam_uy(z_lev) * 1e3, z_lev, "k--", lw=1, label="beam")
plt.plot([beam_uy(zs) * 1e3 for zs in STATIONS.values()], list(STATIONS.values()),
         "o", color="#D55E00", label="stations")
plt.axhspan(1.5, 2.5, color="0.9", label="core")
plt.xlabel(r"u$_y$ (mm)")
plt.ylabel("z (m)")
plt.title("IAEA VP1, constant 15-15Ti: free bow")
plt.grid(True, ls="--", alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(CASE_DIR, "output", "bow_uy.png"), dpi=150)

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)
