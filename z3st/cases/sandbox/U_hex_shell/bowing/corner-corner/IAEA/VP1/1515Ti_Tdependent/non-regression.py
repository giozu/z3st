#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
# Author: Bianca Funaro
"""
Z3ST case: sandbox/U_hex_shell/bowing/corner-corner/IAEA/VP1/1515Ti_Tdependent

IAEA VP1 free thermal bow with the T-dependent 15-15Ti card (eta(T) and E(T)
of the card, Sobolev or DiMarcello, materials/steel_1515Ti_Tdep.py); mesh, BCs
and T as ../ (constant 15-15Ti). The beam of ../non-regression.py is extended
to a thermal strain non-linear in T: in each section z, on the mid-wall hexagon
(corners at ±R on y), eps_th(T(y, z)) is fitted by plane sections,
eps ≈ a + b y, weighted by E(T) ds, and

    kappa(z) = -b(z),   u_y(z) = int_0^z kappa(s) (z - s) ds,

    eps_th = (eta(T) - eta(T_ref)) / (1 + eta(T_ref)),  T_ref = 400 °C.

For constant alpha and E this is ../non-regression.py exactly.
Validity: thin wall (eps_th on the mid-wall line), plane sections.

  - section-mean u_y at the core top, ACLP and duct top against the beam;
  - von Mises maximum and median: no closed form, tracked (gold only).

The gold is compared only under the same hypotheses (card, BCs, geometry,
utils/non_regression.case_hypotheses): a run with another law is NOT COMPARABLE.
"""

import importlib
import os

import matplotlib.pyplot as plt
import numpy as np

from z3st.materials import steel_1515Ti_Tdep as st
from z3st.utils.non_regression import case_hypotheses, finish, load_yaml, metric, tracked
from z3st.utils.utils_extract_xdmf import extract_field_xdmf

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(os.path.abspath(__file__))
inp = load_yaml(CASE_DIR, "input.yaml")
geom = load_yaml(CASE_DIR, inp["geometry_path"])
mat = load_yaml(CASE_DIR, inp["materials"]["steel"])
XDMF_FILE = os.path.join(CASE_DIR, "output", f"{inp['output']['filename']}.xdmf")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")

T_ref = float(mat["T_ref"])        # K
law = mat["expansion_law"]
_module, _func = mat["E"].rsplit(".", 1)
E_card = getattr(importlib.import_module(f"z3st.{_module}"), _func)   # E(T) hook of the card (Pa)
R =(float(geom["D"]) - 2.0 * float(geom["t"]) / np.sqrt(3.0)) / 2.0  # m, mid-wall corner radius

STATIONS = {"core_top": 2.5, "ACLP": 3.0, "duct_top": 4.0}  # m

# z3st sits ~0.015 mm below the beam at every station, 1.4 % of the 1 mm at
# the core top, as in ../non-regression.py.
TOLERANCE = 2.0e-2


# --.. ..- .-.. .-.. --- analytic bow --.. ..- .-.. .-.. ---
def temperature(y, z):
    """T (K) of boundary_conditions.yaml."""
    return 673.15 + np.clip((z - 1.5) / 1.0, 0.0, 1.0) * (125.0 - 25.0 * y / 0.075)


def midwall_hexagon(n_side=400):
    """Points (x, y) and arc-length weights of the mid-wall line, corners on ±y."""
    th = np.radians(90.0 + 60.0 * np.arange(7))
    P = np.stack([R * np.cos(th), R * np.sin(th)], axis=1)
    s = (np.arange(n_side) + 0.5) / n_side
    pts = np.concatenate([P[i] + s[:, None] * (P[i + 1] - P[i]) for i in range(6)])
    return pts[:, 1], np.full(len(pts), R / n_side)


def beam_uy(z):
    y, ds = midwall_hexagon()
    kappa = np.empty_like(z)
    for i, zi in enumerate(z):
        T = temperature(y, zi)
        w = E_card(T) * ds
        e = st.thermal_strain(T, T_ref, law)
        A = [[w.sum(), (w * y).sum()], [(w * y).sum(), (w * y * y).sum()]]
        kappa[i] = -np.linalg.solve(A, [(w * e).sum(), (w * e * y).sum()])[1]
    slope = np.concatenate([[0.0], np.cumsum(0.5 * (kappa[1:] + kappa[:-1]) * np.diff(z))])
    return np.concatenate([[0.0], np.cumsum(0.5 * (slope[1:] + slope[:-1]) * np.diff(z))])


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

z_beam = np.linspace(0.0, float(geom["H"]), 4001)
uy_beam = beam_uy(z_beam)

errors = {}
for key, zs in STATIONS.items():
    num, ref = np.interp(zs, z_lev, uy_mean), np.interp(zs, z_beam, uy_beam)
    print(f"[INFO] u_y at {key:9s} (z = {zs:.1f} m): z3st = {num * 1e3:.4f} mm, beam = {ref * 1e3:.4f} mm")
    errors[f"uy_{key}"] = metric(num, ref)
print(f"[INFO] von Mises: max = {vm.max() * 1e-6:.3f} MPa, median = {np.median(vm) * 1e-6:.3f} MPa")
errors["vonMises_max"] = tracked(vm.max())
errors["vonMises_median"] = tracked(np.median(vm))

# --.. ..- .-.. .-.. --- plot --.. ..- .-.. .-.. ---
plt.figure(figsize=(6, 7))
plt.plot(uy_mean * 1e3, z_lev, "-", color="#0072B2", lw=2, label="z3st, section mean")
plt.plot(uy_beam * 1e3, z_beam, "k--", lw=1, label=r"beam, $\eta(T)$ and $E(T)$")
plt.plot([np.interp(zs, z_beam, uy_beam) * 1e3 for zs in STATIONS.values()], list(STATIONS.values()),
         "o", color="#D55E00", label="stations")
plt.axhspan(1.5, 2.5, color="0.9", label="core")
plt.xlabel(r"u$_y$ (mm)")
plt.ylabel("z (m)")
plt.title("IAEA VP1, T-dependent 15-15Ti: free bow")
plt.grid(True, ls="--", alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(CASE_DIR, "output", "bow_uy.png"), dpi=150)

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
finish(errors, TOLERANCE, OUT_JSON, CASE_DIR, hypotheses=case_hypotheses(CASE_DIR, inp))
