#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
# Author: Bianca Funaro
"""
Z3ST case: sandbox/U_hex_shell/inner_q_outer_T

Hexagonal steel shell (assembly wrapper): heat flux q into the inner faces, T_o
on the outer faces, no mechanical load (rigid-body modes removed). Away from
the corners and the ends each flat is a slab with T = T(xi); the closed section
suppresses the bending, and with k linear in T

    k((T_i + T_o)/2) (T_i - T_o) = q t
    sigma_nn = 0,  sigma_ss = sigma_zz = alpha E / (1 - nu) (T_mean - T(xi))

The comparison is along xi (through the wall), at mid-flat and mid-height,
on the six flats. sigma_zz carries a uniform membrane offset (axial balance
over the whole section, corners included): it is removed before the check
and tracked. Derivation in "Exercise - Thermal stress in hexagone- flats.ipynb".
"""

import importlib
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

from z3st.utils.non_regression import error_metric, finish, load_case, metric, tracked
from z3st.utils.utils_extract_xdmf import extract_field_xdmf

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(os.path.abspath(__file__))
XDMF_FILE = os.path.join(CASE_DIR, "output", "fields.xdmf")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")

geom, inp, mat = load_case(CASE_DIR)
D, t, H = float(geom["D"]), float(geom["t"]), float(geom["H"])  # m
nf, nz = int(geom["nf"]), int(geom["nz"])                       # elements along a flat, along z
E, nu, alpha = float(mat["E"]), float(mat["nu"]), float(mat["alpha"])  # Pa, -, 1/K

sys.path.insert(0, CASE_DIR)
module, func = mat["k"].rsplit(".", 1)
k = getattr(importlib.import_module(module), func)  # W/m·K, k(T) of the card

q = 7000.0   # W/m²  into the inner faces (boundary_conditions.yaml)
To = 500.0   # K     outer faces (boundary_conditions.yaml)

R_o = D / 2                                   # m  outer circumradius
R_i = R_o - 2 * t / np.sqrt(3)                # m  inner circumradius
a_mid = 0.5 * (R_i + R_o) * np.cos(np.pi / 6)  # m  mid-wall apothem
band_s = 0.75 * 0.5 * (R_i + R_o) / nf        # m  half-width of the line, along s
band_z = 0.75 * H / nz                        # m  half-width of the line, along z

TOLERANCE = 1.0e-2  # -


# --.. ..- .-.. .-.. --- analytic functions --.. ..- .-.. .-.. ---
def wall_drop():
    """T_i - T_o from k((T_i + T_o)/2) (T_i - T_o) = q t (exact for linear k)."""
    dT = q * t / k(To)
    for _ in range(50):
        dT = q * t / k(To + dT / 2)
    return dT


def analytic_T(xi):
    """Kirchhoff: int_To^T k dT = q (t/2 - xi), k linear in T."""
    k0, k1 = k(To), k(To + 1.0) - k(To)
    return To + (-k0 + np.sqrt(k0**2 + 2 * k1 * q * (t / 2 - xi))) / k1


def analytic_sigma(xi):
    """sigma_ss = sigma_zz = alpha E / (1 - nu) (T_mean - T), T linear in xi."""
    return alpha * E / (1 - nu) * dT_ref * xi / t


def local_frame(x, y):
    """Outward normal n, tangent s and coordinates xi, s of the flat each point lies on."""
    phi = np.floor(np.arctan2(y, x) / (np.pi / 3)) * np.pi / 3 + np.pi / 6
    n = np.stack([np.cos(phi), np.sin(phi), 0 * phi], axis=1)
    s = np.stack([-np.sin(phi), np.cos(phi), 0 * phi], axis=1)
    p = np.stack([x, y, 0 * x], axis=1)
    return n, s, (p * n).sum(1) - a_mid, (p * s).sum(1)


def through_wall(xi, values):
    """Mean of ``values`` per xi layer: the six flats collapse onto one curve."""
    layers, inv = np.unique(np.round(xi, 9), return_inverse=True)
    return layers, np.bincount(inv, weights=values) / np.bincount(inv)


def rel_l2(num, ref, scale=None):
    scale = np.sqrt(np.mean(ref**2)) if scale is None else scale
    return float(np.sqrt(np.mean((num - ref) ** 2)) / scale)


# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
dT_ref = wall_drop()

x, y, z, T_all = extract_field_xdmf(XDMF_FILE, field_name="Temperature", step_index=-1)
_, _, xi, s = local_frame(x, y)
line = (np.abs(s) < band_s) & (np.abs(z - H / 2) < band_z)
xi_T, T = through_wall(xi[line], T_all[line])
T_ref = analytic_T(xi_T)

xc, yc, zc, S = extract_field_xdmf(XDMF_FILE, field_name="Stress", step_index=-1)
S = S.reshape(-1, 3, 3)
n, tg, xi_c, s_c = local_frame(xc, yc)
line = (np.abs(s_c) < band_s) & (np.abs(zc - H / 2) < band_z)
xi_s, s_ss = through_wall(xi_c[line], np.einsum("ci,cij,cj->c", tg, S, tg)[line])
_, s_zz = through_wall(xi_c[line], S[line, 2, 2])
_, s_nn = through_wall(xi_c[line], np.einsum("ci,cij,cj->c", n, S, n)[line])
sigma_ref = analytic_sigma(xi_s)
zz_offset = s_zz.mean()

print(f"[INFO] Wall drop T_i - T_o : z3st = {T[0] - T[-1]:.5f} K, analytic = {dT_ref:.5f} K")
for j in range(len(xi_s)):
    print(f"[INFO] xi = {xi_s[j] * 1e3:+.3f} mm : analytic = {sigma_ref[j] * 1e-6:+.4f}, "
          f"ss = {s_ss[j] * 1e-6:+.4f}, zz = {s_zz[j] * 1e-6:+.4f}, nn = {s_nn[j] * 1e-6:+.4f} MPa")
print(f"[INFO] sigma_zz membrane offset : {zz_offset * 1e-6:+.4f} MPa")

# --.. ..- .-.. .-.. --- plot --.. ..- .-.. .-.. ---
xi_a = np.linspace(-t / 2, t / 2, 101)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
ax1.plot(xi_T * 1e3, T, "o", color="#0072B2", label="z3st")
ax1.plot(xi_a * 1e3, analytic_T(xi_a), "k--", label="analytic")
ax1.set_xlabel(r"$\xi$ (mm)")
ax1.set_ylabel("T (K)")
for vals, col, lab in [(s_ss, "#009E73", r"$\sigma_{ss}$"), (s_zz - zz_offset, "#D55E00", r"$\sigma_{zz}$ - offset"),
                       (s_nn, "#CC79A7", r"$\sigma_{nn}$")]:
    ax2.plot(xi_s * 1e3, vals * 1e-6, "o", color=col, label=lab)
ax2.plot(xi_a * 1e3, analytic_sigma(xi_a) * 1e-6, "k--", label=r"analytic $\sigma_{ss} = \sigma_{zz}$")
ax2.set_xlabel(r"$\xi$ (mm)")
ax2.set_ylabel("Stress (MPa)")
for ax in (ax1, ax2):
    ax.grid(True, ls="--", alpha=0.6)
    ax.legend()
fig.suptitle(f"Through the wall, mid-flat and mid-height: q = {q:.0f} W/m², T_o = {To:.0f} K")
fig.tight_layout()
fig.savefig(os.path.join(CASE_DIR, "output", "through_wall.png"), dpi=150)

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
L2_T = float(np.sqrt(np.mean((T - T_ref) ** 2)))
Linf_T = float(np.max(np.abs(T - T_ref)))
scale = np.sqrt(np.mean(sigma_ref**2))

errors = {
    "L2_error_T": error_metric(L2_T, rel=L2_T / dT_ref),
    "Linf_error_T": error_metric(Linf_T, rel=Linf_T / dT_ref),
    "T_max": metric(T.max(), T_ref.max()),
    "wall_dT": metric(T[0] - T[-1], dT_ref),
    "L2_error_sigma_ss": error_metric(rel_l2(s_ss, sigma_ref)),
    "L2_error_sigma_zz": error_metric(rel_l2(s_zz - zz_offset, sigma_ref)),
    "L2_error_sigma_nn": error_metric(rel_l2(s_nn, 0 * s_nn, scale)),
    "sigma_zz_offset": tracked(zz_offset),
}

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)
