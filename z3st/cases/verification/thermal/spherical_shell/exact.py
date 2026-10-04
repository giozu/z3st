#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""Exact solution of the spherical_shell case, shared by non-regression.py and
plots.py.

Steady conduction in a spherical shell Ri <= r <= Ro with the gamma source of
core/spine.py for geometry_type sphere,

    q(r) = q0 (Rg / r) exp(-mu (r - Rg)),     Rg = gamma_inner_radius or Ri,

and fixed temperatures Ti (inner) and To (outer). The heat equation
(1/r^2) d/dr(r^2 k dT/dr) = -q integrates in closed form, because
d/dr[exp(-mu (r - Rg)) / r] supplies both terms of r^2 dT/dr:

    T(r) = C2 - C1/r - (q0 Rg / (k mu^2)) exp(-mu (r - Rg)) / r.

The thermo-elastic state for that T(r), with a free inner surface and a
clamped outer surface, is

    u(r) = (1+nu)/(1-nu) alpha/r^2 int_Ri^r (T - T_sf) s^2 ds + C3 r + C4/r^2,

with sigma_rr(Ri) = 0 and u(Ro) = 0. T_sf is the stress-free temperature.

Every parameter is read from the case files, so the reference cannot drift
from the input.
"""

import os

import numpy as np
import yaml

from z3st.utils.non_regression import load_case


def parameters(case_dir):
    """Geometry, material and boundary values of the case, as a dict."""
    geom, _, mat = load_case(case_dir)
    with open(os.path.join(case_dir, "boundary_conditions.yaml")) as f:
        bc = {b["region"]: float(b["temperature"])
              for b in yaml.safe_load(f)["thermal"]["mat0"]}
    Ri = float(geom["Ri"])
    return {
        "Ri": Ri,
        "Ro": float(geom["Ro"]),
        "k": float(mat["k"]),
        "q0": float(mat["gamma_heating"]),
        "mu": float(mat["mu_gamma"]),
        "Rg": float(mat.get("gamma_inner_radius", Ri)),
        "Ti": bc["inner"],
        "To": bc["outer"],
        "E": float(mat["E"]),
        "nu": float(mat["nu"]),
        "alpha": float(mat["alpha"]),
        "T_sf": float(mat["T_ref"]),
    }


def source(r, p):
    """Gamma heating q(r) (W/m3), as implemented in core/spine.py."""
    return p["q0"] * (p["Rg"] / r) * np.exp(-p["mu"] * (r - p["Rg"]))


def temperature(r, p):
    """Exact T(r) of the spherical shell (K)."""
    def g(rr):
        return -(p["q0"] * p["Rg"] / (p["k"] * p["mu"] ** 2)) * np.exp(-p["mu"] * (rr - p["Rg"])) / rr

    # T = C2 - C1/r + g(r): two linear conditions for C1, C2
    A = np.array([[-1.0 / p["Ri"], 1.0], [-1.0 / p["Ro"], 1.0]])
    C1, C2 = np.linalg.solve(A, np.array([p["Ti"] - g(p["Ri"]), p["To"] - g(p["Ro"])]))
    return C2 - C1 / r + g(r)


def temperature_slab(r, p):
    """Plane-slab approximation the case used before: exponential source
    q0 exp(-mu x), x = r - Ri, in a slab of thickness Ro - Ri. Kept only to
    show, in plots.py, how far it is from the sphere."""
    L = p["Ro"] - p["Ri"]
    x = r - p["Ri"]
    mu, k, q0 = p["mu"], p["k"], p["q0"]
    return (p["Ti"] + (p["To"] - p["Ti"]) * x / L
            + q0 / (mu**2 * k) * ((x / L) * (np.exp(-mu * L) - 1.0) - (np.exp(-mu * x) - 1.0)))


def mechanics(r_eval, p, n_grid=20001):
    """Exact u_r (m), sigma_rr and sigma_tt (Pa) at the radii r_eval."""
    E, nu, alpha = p["E"], p["nu"], p["alpha"]
    rg = np.linspace(p["Ri"], p["Ro"], n_grid)
    dT = temperature(rg, p) - p["T_sf"]
    f = dT * rg**2
    I = np.concatenate(([0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(rg))))
    kk = (1.0 + nu) / (1.0 - nu) * alpha
    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    G = E / (2.0 * (1.0 + nu))
    bulk3 = E / (1.0 - 2.0 * nu)  # 3 lambda + 2 G

    def fields(C3, C4):
        u = kk * I / rg**2 + C3 * rg + C4 / rg**2
        du = kk * (dT - 2.0 * I / rg**3) + C3 - 2.0 * C4 / rg**3
        tr = du + 2.0 * u / rg
        srr = lam * tr + 2.0 * G * du - bulk3 * alpha * dT
        stt = lam * tr + 2.0 * G * u / rg - bulk3 * alpha * dT
        return u, srr, stt

    # sigma_rr(Ri) and u(Ro) are linear in (C3, C4)
    base = fields(0.0, 0.0)
    e3 = [a - b for a, b in zip(fields(1.0, 0.0), base)]
    e4 = [a - b for a, b in zip(fields(0.0, 1.0), base)]
    A = np.array([[e3[1][0], e4[1][0]], [e3[0][-1], e4[0][-1]]])
    C3, C4 = np.linalg.solve(A, -np.array([base[1][0], base[0][-1]]))
    u, srr, stt = fields(C3, C4)
    return (np.interp(r_eval, rg, u), np.interp(r_eval, rg, srr), np.interp(r_eval, rg, stt))


def demo():
    """Self-check: the exact T satisfies the boundary values and the heat
    equation, and the mechanics satisfies its two boundary conditions."""
    p = {"Ri": 2.0, "Ro": 2.5, "k": 48.1, "q0": 2.0e6, "mu": 24.0, "Rg": 2.0,
         "Ti": 494.0, "To": 491.0, "E": 1.77e11, "nu": 0.3, "alpha": 1.7e-5, "T_sf": 300.0}
    assert abs(temperature(np.array(2.0), p) - 494.0) < 1e-9
    assert abs(temperature(np.array(2.5), p) - 491.0) < 1e-9
    r = np.linspace(2.05, 2.45, 9)
    h = 1e-3
    # (1/r^2) d/dr(r^2 dT/dr) by central differences on the half-step radii
    lap = ((r + h / 2) ** 2 * (temperature(r + h, p) - temperature(r, p))
           - (r - h / 2) ** 2 * (temperature(r, p) - temperature(r - h, p))) / (h**2 * r**2)
    assert np.allclose(p["k"] * lap + source(r, p), 0.0, atol=1e-3 * source(r, p).max())
    u, srr, _ = mechanics(np.array([2.0, 2.5]), p)
    assert abs(srr[0]) < 1e-6 * p["E"] * 1e-3 and abs(u[1]) < 1e-15
    print("exact.py: self-check passed")


if __name__ == "__main__":
    demo()
