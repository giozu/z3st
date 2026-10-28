#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Constants and analytic slab solution of this case, shared by ``plots.py`` and
the notebook. Everything is read from the YAML files named in ``input.yaml``;
the only default is ``orientation: flat``, as in ``mesh.geo``.
``non-regression.py`` does not import it (self-contained, as the native cases).

Local frame of a flat: xi through the wall (0 at mid-wall, outward), s along
the flat (0 at mid-flat), z along the axis.
"""

import importlib
import os
import sys

import numpy as np
import yaml

from z3st.utils.utils_extract_xdmf import extract_field_xdmf

HERE = os.path.dirname(os.path.abspath(__file__))
XDMF = os.path.join(HERE, "output", "fields.xdmf")


def _load(name):
    with open(os.path.join(HERE, name), "r") as fh:
        return yaml.safe_load(fh)


_input = _load("input.yaml")
_geom = _load(_input["geometry_path"])
_bcs = _load(_input["boundary_conditions_path"])
_mat = _load(_input["materials"]["steel"])

# --- Geometry (m) -----------------------------------------------------------
D = float(_geom["D"])                               # outer circumscribed diameter
T_WALL = float(_geom["t"])                          # wall thickness
H = float(_geom["H"])                               # height
ORIENTATION = _geom.get("orientation", "flat")
R_O = D / 2                                         # outer circumradius
R_I = R_O - 2 * T_WALL / np.sqrt(3)                 # inner circumradius
A_MID = 0.5 * (R_I + R_O) * np.cos(np.pi / 6)       # mid-wall apothem
L_MID = 0.5 * (R_I + R_O)                           # flat length at mid-wall
A0 = np.pi / 2 if ORIENTATION == "vertical" else 0.0

# --- Mesh (transfinite, see mesh.geo) ---------------------------------------
NF, NT, NZ = int(_geom["nf"]), int(_geom["nt"]), int(_geom["nz"])
H_S = L_MID / NF                                    # element width at mid-wall
H_Z = H / NZ                                        # element height
XI_CELLS = -T_WALL / 2 + (np.arange(NT) + 0.5) * T_WALL / NT  # cell-centre layers
BAND = 0.75                                         # sampling half-width, in elements

# --- Material ---------------------------------------------------------------
E = float(_mat["E"])                                # (Pa)
NU = float(_mat["nu"])                              # (-)
ALPHA = float(_mat["alpha"])                        # (1/K)
T_REF = float(_mat["T_ref"])                        # (K) stress-free temperature
_module, _func = _mat["k"].split(".")
sys.path.insert(0, HERE)
k = getattr(importlib.import_module(_module), _func)  # (W/m.K)

# --- Boundary conditions ----------------------------------------------------
_thermal = {bc["type"]: bc for bc in _bcs["thermal"]["steel"]}
Q = -float(_thermal["Neumann"]["flux"])             # (W/m^2) into the wall
T_O = float(_thermal["Dirichlet"]["temperature"])   # (K) outer face


# --- Analytic slab solution -------------------------------------------------
def wall_drop():
    """T_i - T_o from k((T_i + T_o)/2) (T_i - T_o) = q t, exact for linear k."""
    dT = Q * T_WALL / k(T_O)
    for _ in range(50):
        dT = Q * T_WALL / k(T_O + dT / 2)
    return dT


DT = wall_drop()                                    # (K)
T_I = T_O + DT                                      # (K) inner face
C_TH = ALPHA * E / (1 - NU)                         # (Pa/K) biaxial thermal modulus
SIGMA_SURF = C_TH * DT / 2                          # (Pa) surface stress


def temperature(xi):
    """Kirchhoff solution: integral of k from T_o to T(xi) equals q (t/2 - xi)."""
    k_o = k(T_O)
    k_1 = k(T_O + 1.0) - k_o
    return T_O + (-k_o + np.sqrt(k_o**2 + 2 * k_1 * Q * (T_WALL / 2 - xi))) / k_1


def sigma_wall(xi):
    """sigma_ss = sigma_zz = alpha E / (1 - nu) (T_mean - T), for linear T(xi)."""
    return C_TH * DT * np.asarray(xi) / T_WALL


# --- z3st fields in the local frame -----------------------------------------
def flat_frame(x, y):
    """Normal, tangent, xi and s of each point, from the flat it lies on."""
    sector = np.floor((np.arctan2(y, x) - A0) / (np.pi / 3))
    phi_n = A0 + sector * np.pi / 3 + np.pi / 6
    zero = np.zeros_like(phi_n)
    n = np.stack([np.cos(phi_n), np.sin(phi_n), zero], axis=1)
    tang = np.stack([-np.sin(phi_n), np.cos(phi_n), zero], axis=1)
    p = np.stack([x, y, np.zeros_like(x)], axis=1)
    return n, tang, (p * n).sum(1) - A_MID, (p * tang).sum(1)


def local_fields(xdmf=XDMF, step_index=-1):
    """Nodal temperature and cell stress (DG0) of a step, in the local frame.

    Returns two dicts: ``nodes`` (xi, s, z, T) and ``cells``
    (xi, s, z, nn, ss, zz).
    """
    x, y, z, T = extract_field_xdmf(xdmf, "Temperature", step_index=step_index)
    _, _, xi, s = flat_frame(x, y)
    nodes = {"xi": xi, "s": s, "z": z, "T": T}

    xc, yc, zc, S = extract_field_xdmf(xdmf, "Stress", step_index=step_index)
    S = S.reshape(-1, 3, 3)
    n, tang, xi_c, s_c = flat_frame(xc, yc)
    cells = {
        "xi": xi_c, "s": s_c, "z": zc,
        "nn": np.einsum("ci,cij,cj->c", n, S, n),
        "ss": np.einsum("ci,cij,cj->c", tang, S, tang),
        "zz": S[:, 2, 2],
    }
    return nodes, cells


def through_wall(f):
    """Line along xi: mid-flat, mid-height."""
    return (np.abs(f["s"]) < BAND * H_S) & (np.abs(f["z"] - H / 2) < BAND * H_Z)


def along_flat(f):
    """Line along s: mid-height, every flat."""
    return np.abs(f["z"] - H / 2) < BAND * H_Z


def along_axis(f):
    """Line along z: mid-flat, every flat."""
    return np.abs(f["s"]) < BAND * H_S


def profile(f, coord, field, mask):
    """Mean of ``field`` over the points of ``mask`` sharing a ``coord`` value.

    The six flats and the one or two rows inside the band collapse onto one
    curve, sorted by ``coord``.
    """
    c = np.round(f[coord][mask], 9)
    axis, inverse = np.unique(c, return_inverse=True)
    values = np.bincount(inverse, weights=f[field][mask]) / np.bincount(inverse)
    return axis, values


def on_layer(f, xi):
    """Points (nodes or cell centres) on the wall layer ``xi``."""
    return np.abs(f["xi"] - xi) < 0.1 * T_WALL / NT


def check_consistency(cells=None):
    """Fail loudly on the inconsistencies that silently invalidate the overlay.

    Returns the list of problems found (empty when the case is coherent).
    """
    problems = []

    if ORIENTATION not in ("flat", "vertical"):
        problems.append(f"orientation is '{ORIENTATION}', expected 'flat' or 'vertical'")

    # The exact wall drop and temperature() assume k linear in T.
    T_mid = T_O + DT / 2
    curvature = k(T_mid + 1.0) - 2 * k(T_mid) + k(T_mid - 1.0)
    if abs(curvature) > 1e-9 * abs(k(T_mid)):
        problems.append("k(T) is not linear: the Kirchhoff closed form is not exact")

    if _input.get("regime") != "3d":
        problems.append(f"regime is '{_input.get('regime')}', the analysis needs '3d'")

    # Every cell centre lies on one of the NT layers when geometry.yaml
    # (orientation, D, t, nt) matches the solved mesh.
    if cells is not None:
        offset = np.abs(cells["xi"][:, None] - XI_CELLS).min(axis=1).max()
        if offset > 0.1 * T_WALL / NT:
            problems.append(
                "cell centres are not on the nt wall layers: geometry.yaml "
                "(orientation, D, t, nt) does not match the solved mesh"
            )

    return problems


def report():
    """Print the resolved parameter set and the analytic reference values."""
    print("[case_params] resolved from the case YAML files:")
    for key in ("D", "T_WALL", "H", "ORIENTATION", "R_O", "R_I", "NF", "NT", "NZ",
                "E", "NU", "ALPHA", "T_REF", "Q", "T_O"):
        print(f"  {key:<12} = {globals()[key]}")
    print(f"  {'DT':<12} = {DT:.6f} K (T_i - T_o)")
    print(f"  {'SIGMA_SURF':<12} = {SIGMA_SURF / 1e6:.4f} MPa")
