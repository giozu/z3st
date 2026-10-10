#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Single source of truth for the post-processing of this case.

Why this module exists: four consumers need the same constants and the same
analytic solution, and each reads them at a different moment:

  * ``diagnostics.py``   inside the z3st run, after every converged step;
  * ``non-regression.py`` after the run, from ``output/history.csv``;
  * ``plots.py``          after the run, from ``output/fields_<regime>.xdmf``;
  * the notebook          interactively.

Defining the geometry, the sample stations and the analytic solution once here
keeps the four from drifting apart: a change to D, t, the mesh counts, the
material card or the boundary conditions reaches all of them through the YAML
files named in ``input.yaml``. The only default is ``orientation: flat``, the
same as in ``mesh.geo``; any other missing key fails here instead of silently
using a value the mesh or the solver did not use.

Two variants, chosen by ``./Allrun 2d|3d`` (``regime`` in ``input.yaml``):
``3d`` is the shell of height H; ``2d`` is its cross-section in plane strain,
with no z line (H = 0, one element "along z").

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


def _load(name):
    with open(os.path.join(HERE, name), "r") as fh:
        return yaml.safe_load(fh)


_input = _load("input.yaml")
# output.filename in input.yaml (fields_2d | fields_3d), same rule as the z3st writer
_fields = os.path.splitext(_input.get("output", {}).get("filename", "fields"))[0]
XDMF = os.path.join(HERE, "output", _fields + ".xdmf")
_geom = _load(_input["geometry_path"])
_bcs = _load(_input["boundary_conditions_path"])
_mat = _load(_input["materials"]["steel"])

# --- Geometry (m) -----------------------------------------------------------
D = float(_geom["D"])                               # outer circumscribed diameter
T_WALL = float(_geom["t"])                          # wall thickness
REGIME = _input["regime"]                           # "2d" | "3d"
IS_3D = REGIME == "3d"
H = float(_geom["H"]) if IS_3D else 0.0             # height
ORIENTATION = _geom.get("orientation", "flat")
R_O = D / 2                                         # outer circumradius
R_I = R_O - 2 * T_WALL / np.sqrt(3)                 # inner circumradius
A_MID = 0.5 * (R_I + R_O) * np.cos(np.pi / 6)       # mid-wall apothem
L_MID = 0.5 * (R_I + R_O)                           # flat length at mid-wall
A0 = np.pi / 2 if ORIENTATION == "vertical" else 0.0

# --- Mesh (transfinite, see mesh.geo) ---------------------------------------
NF, NT = int(_geom["nf"]), int(_geom["nt"])
NZ = int(_geom["nz"]) if IS_3D else 1
H_S = L_MID / NF                                    # element width at mid-wall
H_Z = H / NZ if IS_3D else np.inf                   # element height (2d: no z band)
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


# sigma_zz = sigma_wall + a uniform offset. In 3d the offset comes from the
# axial force balance over the whole section, corners included: no closed form.
# In 2d (plane strain, eps_zz = 0) sigma_zz = nu sigma_ss - E alpha (T - T_ref),
# and with sigma_ss = alpha E / (1 - nu) (T_mean - T) the offset is exact.
SIGMA_ZZ_OFFSET_2D = -E * ALPHA * (T_O + DT / 2 - T_REF)   # (Pa)


# --- Sample stations (read every step by diagnostics.py) ---------------------
def flat_axes(m):
    """Outward normal and tangent of flat m (0..5)."""
    phi_n = A0 + np.pi / 6 + m * np.pi / 3
    return (np.array([np.cos(phi_n), np.sin(phi_n), 0.0]),
            np.array([-np.sin(phi_n), np.cos(phi_n), 0.0]))


def flat_point(m, xi, u_s, u_z):
    """Point of flat m at wall position xi, along-flat fraction u_s (0 and 1
    are the corners) and height fraction u_z."""
    n, tang = flat_axes(m)
    a = A_MID + xi
    s = (u_s - 0.5) * 2 * a * np.tan(np.pi / 6)  # flat length at apothem a
    return a * n + s * tang + np.array([0.0, 0.0, u_z * H])


def _centre(u, n):
    """Fraction of the centre of the cell containing u in [0, 1], n cells."""
    return (min(int(np.floor(u * n)), n - 1) + 0.5) / n


# Boundary points are pulled inside the wall by this fraction of t, so the
# bounding-box search finds them despite round-off.
_EPS = 1e-6
XI_I, XI_O = -T_WALL / 2 * (1 - _EPS), T_WALL / 2 * (1 - _EPS)

# A station is a fixed probe position, given once in the local frame and
# placed on all six flats; diagnostics.py evaluates it after every step and
# writes the flat average to history.csv, where non-regression.py reads it.
# Unrelated to ``thermal.analysis: stationary``.
#
# name: (kind, xi, u_s, u_z) -> history.csv column(s)
#   kind "node": temperature, read at the point          -> <name>_K
#   kind "cell": DG0 stress, at a cell centre (xi on       -> ss_/zz_/nn_<name>_MPa
#                XI_CELLS, u from _centre) so the read
#                cannot fall on a cell boundary
#
#   T_i, T_o              faces, mid-flat, mid-height: wall drop
#   T_i_corner/_bottom/_top  inner face at corner/ends: T uniformity along s, z
#   L0 .. L{nt-1}         every wall layer, mid-flat, mid-height: sigma(xi)
#   window_s, window_z    outer layer at L/4, H/4: edge of the checked far field
#   corner, bottom, top   outer layer at corner/ends: disturbances, tracked only
#
# The window stations sit at the cell containing the point a quarter of the
# line from one end; with u = 0.25 the cell centre lies just inside the
# checked far field (|s| < L/4, |z - H/2| < H/4). By the mirror symmetry
# about mid-flat and mid-height the side does not matter.
STATIONS = {
    "T_i": ("node", XI_I, 0.5, 0.5),
    "T_o": ("node", XI_O, 0.5, 0.5),
    "T_i_corner": ("node", XI_I, _EPS, 0.5),
    **{f"L{j}": ("cell", xi, _centre(0.5, NF), _centre(0.5, NZ))
       for j, xi in enumerate(XI_CELLS)},
    "window_s": ("cell", XI_CELLS[-1], _centre(0.25, NF), _centre(0.5, NZ)),
    "corner": ("cell", XI_CELLS[-1], _centre(1.0, NF), _centre(0.5, NZ)),
}
if IS_3D:  # the z line
    STATIONS.update({
        "T_i_bottom": ("node", XI_I, 0.5, _EPS),
        "T_i_top": ("node", XI_I, 0.5, 1 - _EPS),
        "window_z": ("cell", XI_CELLS[-1], _centre(0.5, NF), _centre(0.25, NZ)),
        "bottom": ("cell", XI_CELLS[-1], _centre(0.5, NF), _centre(0.0, NZ)),
        "top": ("cell", XI_CELLS[-1], _centre(0.5, NF), _centre(1.0, NZ)),
    })


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

    if REGIME not in ("2d", "3d"):
        problems.append(f"regime is '{REGIME}', the analysis needs '2d' or '3d'")

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
    for key in ("REGIME", "D", "T_WALL", "H", "ORIENTATION", "R_O", "R_I", "NF", "NT", "NZ",
                "E", "NU", "ALPHA", "T_REF", "Q", "T_O"):
        print(f"  {key:<12} = {globals()[key]}")
    print(f"  {'DT':<12} = {DT:.6f} K (T_i - T_o)")
    print(f"  {'SIGMA_SURF':<12} = {SIGMA_SURF / 1e6:.4f} MPa")
