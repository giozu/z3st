#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.4.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""Independent radial reference for the relaxation of creep_shrink_fit_2D.

The same joint as the case, reduced to the radius: an elastic solid shaft
(the pellet, axially free, sigma_r = sigma_t = -P), a hub (the cladding)
creeping by Norton's law with the J2 flow rule, axially free in generalised
plane strain (uniform eps_z, zero axial force), and the penalty spring of the
contact model in series at the interface. The hub is discretised by linear
elements in r with one integration point, and the creep strain is advanced
explicitly with a step that keeps each increment below a fraction of the
elastic strain.

It shares no code with Z3ST and no approximation with Esposito eq. (21): no
Tresca criterion and no assumed kinematics of the creep strain rate. It reads
every constant through case_params.

Converged values (dt_frac = 5e-4, 200 elements): 24.700, 7.120, 5.063 and
3.605 MPa at 0, 600, 1240 and 2500 days. Halving dt_frac moves them by less
than 0.1 %.

Run: python3 reference_1d.py   (self-check)
"""

import numpy as np

import case_params as cp


def pressure_history(days, n_el=200, dt_frac=5.0e-4, T=580.0):
    """Contact pressure (Pa) at the requested times (days)."""
    b1, b2, c = cp.R_PELLET, cp.R_CLAD_I, cp.R_CLAD_O
    E1, nu1, E2, nu2 = cp.E_FUEL, cp.NU_FUEL, cp.E_CLAD, cp.NU_CLAD
    A = cp.CREEP_A0 * np.exp(-cp.CREEP_Q / (cp.R_GAS * T))
    n = cp.CREEP_N
    delta = cp.hot_interference(T)
    # shaft compliance (solid, axially free) and penalty spring in series
    keff = 1.0 / (1.0 / cp.K_PEN + b1 * (1.0 - nu1) / E1)

    r = np.linspace(b2, c, n_el + 1)
    rc, h = 0.5 * (r[1:] + r[:-1]), np.diff(r)
    lam = E2 * nu2 / ((1 + nu2) * (1 - 2 * nu2))
    G = E2 / (2 * (1 + nu2))
    C = lam * np.ones((3, 3)) + 2 * G * np.eye(3)
    ndof = n_el + 2                    # nodal u_r and one uniform eps_z

    # strain at the element centre: [eps_r, eps_t, eps_z] = B u
    B = np.zeros((n_el, 3, ndof))
    e = np.arange(n_el)
    B[e, 0, e], B[e, 0, e + 1] = -1.0 / h, 1.0 / h
    B[e, 1, e], B[e, 1, e + 1] = 0.5 / rc, 0.5 / rc
    B[:, 2, -1] = 1.0
    w = rc * h
    K = np.einsum("eia,ij,ejb,e->ab", B, C, B, w)
    K[0, 0] += keff * b2
    K_inv = np.linalg.inv(K)

    eps_cr = np.zeros((n_el, 3))
    t, t_end = 0.0, max(days) * 86400.0
    targets = sorted(d * 86400.0 for d in days)
    out = {}
    dt = 1.0
    while True:
        F = np.einsum("eia,ij,ej,e->a", B, C, eps_cr, w)
        F[0] += keff * delta * b2
        u = K_inv @ F
        P = keff * (delta - u[0])
        while targets and t >= targets[0] - 1e-6:
            out[targets.pop(0)] = P
        if not targets:
            break
        sig = np.einsum("ij,ej->ei", C, np.einsum("eia,a->ei", B, u) - eps_cr)
        s = sig - sig.mean(axis=1, keepdims=True)
        seq = np.sqrt(1.5 * (s**2).sum(axis=1))
        dt_max = dt_frac * (seq.max() / E2) / (A * seq.max() ** n)
        dt = min(dt_max, 1.2 * dt, targets[0] - t)
        eps_cr += dt * 1.5 * A * seq[:, None] ** (n - 1) * s
        t += dt
    return np.array([out[d * 86400.0] for d in days])


def demo():
    """Self-check: at t = 0 the reference is the elastic Lame joint with the
    penalty spring in series, which case_params.elastic_factor gives in closed
    form."""
    p0 = pressure_history([0.0])[0]
    f = cp.elastic_factor(0.0, cp.R_CLAD_I, cp.R_CLAD_O, cp.E_FUEL, cp.NU_FUEL, cp.E_CLAD, cp.NU_CLAD)
    p_lame = cp.hot_interference(580.0) / (1.0 / f + 1.0 / cp.K_PEN)
    assert abs(p0 - p_lame) / p_lame < 2e-3, (p0, p_lame)
    days = [0.0, 600.0, 1240.0, 2500.0]
    p = pressure_history(days)
    assert np.all(np.diff(p) < 0.0)
    print("reference_1d: self-check passed")
    for d, pi in zip(days, p):
        print(f"  {d:6.0f} d   {pi / 1e6:.3f} MPa")


if __name__ == "__main__":
    demo()
