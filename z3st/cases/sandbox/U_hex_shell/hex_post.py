#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Common post-processing of the U_hex_shell subfolders (2d, 3d, IAEA_VP1,
T_asymmetries, ...), read only from a z3st XDMF/HDF5 output.

Usage::

    python hex_post.py <subfolder> [<subfolder> ...]

Each subfolder carries a ``post.yaml`` (paths relative to the subfolder)::

    fields: output/fields_3d.xdmf      # z3st output
    material: ../15_15Ti_local.yaml    # E, nu, alpha, T_ref (and k for 'radial'); 'vp1' also E(T) + eigenstrain hooks
    temperature: field                 # field | module.func, f(x) with x = (x, y, z)
    support: free                      # clamped (bottom fixed) | free
    closed_form: none                  # none | radial | vp1:<module> (<module> may be a path, ../vp1_temperature)
    boundary_conditions: ../boundary_conditions.yaml   # 'radial' only: q and T_o
    stations_z: [0.25, 0.5, 0.75, 1.0] # heights of the tables (m), optional
    along_s_z: {inlet core: 1.5, tip: 4.0}  # heights of the along-s figures (m), optional

Every figure goes to ``results/`` next to this file, named
``<subfolder>_<figure>.png`` and titled with the subfolder name. The first
two follow the IAEA VP1 page of the INL Virtual Test Bed (hex_duct_linear):

  vtb_temperature_bowing  3D view: temperature, and the deformed duct
                          (magnified) coloured by the displacement along the
                          bowing direction
  vtb_bowing              bowing displacement (x) against the axial position
                          (y): FEM as lines, closed form as circles, and a
                          table at the stations
  curvature               curvature against the axial position, same layout
  displacement_z          section-mean displacement along z
  stress_z                max von Mises and max |sigma_ss| of each section
  along_s_<tag>           von Mises, hoop sigma_ss, strains and u_n along the
                          perimeter, with the vertices marked and the closed
                          forms (equivalent cylinder, hexagon) where they exist,
                          at the max von Mises and at mid-height
  along_s_<group>         with ``along_s_z`` instead: stress (sigma_ss, von
                          Mises), strain_sz (eps_ss, eps_zz), normal (eps_nn,
                          u_n), one row per height and a recap row with the
                          heights overlaid
  max_stress              table: magnitude and position of the max stresses

Closed forms (``closed_form``), never fitted to the FEM:

  radial  through-wall gradient (inner flux q, outer T_o, linear k(T)).
          No bow by the D6 symmetry (kappa = w = 0, shown as a table).
          Along s: thick cylinder of the same perimeter and thickness
          (Timoshenko; plane strain in 2d, free ends in 3d, far from the
          ends), and the slab sigma_ss in the middle half of each flat.
  vp1     T linear across the section, T1 / T2 the cold / hot corner
          temperatures of <module>: kappa = alpha (T2 - T1) / D towards the
          cold corner, w from u = u' = 0 at the clamp, u_z from the section
          mean temperature. Where T is linear in (x, y, z) the free solution
          is stress-free with u = alpha [theta x - |x|^2 grad T / 2], valid
          away from the clamp, for the hexagon and the cylinder alike.
          With a T-dependent card (E(T) hook, eigenstrain hook in place of
          alpha): kappa and u_z from the plane-section fit of eps_th(T) on
          the mid-wall line, weighted by E(T) ds; FEM only along s.
  none    FEM only, said so in the figure.

With several runs, each run is compared with its own closed form.

Local frame of a flat: xi through the wall (0 at mid-wall, outward), s along
the flat. The perimeter coordinate p runs at mid-wall from the first vertex,
counter-clockwise, flat after flat.
"""

import importlib
import json
import os
import sys

import h5py
import matplotlib
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import yaml  # noqa: E402

os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")   # results and comparisons of every subfolder
COMPASS = {90: "nord", 30: "nord-est", 330: "sud-est", 270: "sud", 210: "sud-west", 150: "nord-west"}
RUN_COLORS = ("tab:blue", "tab:orange", "tab:green", "tab:purple")
LAYER_STYLES = {"inner": "-", "outer": "--"}


# --. reading --..
def _last(group, step_index):
    keys = sorted(group, key=lambda k: float(k.replace("_", ".")))
    return group[keys[step_index]][:]


def _tensor(a):
    """(n, k*k) -> (n, 3, 3), a 2x2 tensor padded with zeros."""
    k = int(round(np.sqrt(a.shape[1])))
    out = np.zeros((len(a), 3, 3))
    out[:, :k, :k] = a.reshape(-1, k, k)
    return out


def _pad3(a):
    out = np.zeros((len(a), 3))
    out[:, :a.shape[1]] = a
    return out


def _import(name, case_dir):
    """Module ``name`` of ``case_dir``: two case folders may hold a module of the
    same name (VP1/vp1_temperature.py and VP1/z_uniform/vp1_temperature.py)."""
    case_dir = os.path.abspath(case_dir)
    if case_dir in sys.path:
        sys.path.remove(case_dir)
    sys.path.insert(0, case_dir)
    mod = sys.modules.get(name)
    if mod is not None and os.path.dirname(os.path.abspath(getattr(mod, "__file__", "") or "")) != case_dir:
        del sys.modules[name]
    return importlib.import_module(name)


def load_z3st(xdmf, step_index=-1, temperature="field", case_dir=HERE):
    """Mesh and fields of a step: X, cells (VTK order), u (nodes), T (nodes),
    S and E (cell 3x3 tensors, DG0)."""
    with h5py.File(os.path.splitext(xdmf)[0] + ".h5", "r") as h:
        X = _pad3(h["Mesh/mesh/geometry"][:])
        cells = h["Mesh/mesh/topology"][:]
        f = h["Function"]
        u = _pad3(_last(f["Displacement"], step_index))
        S = _tensor(_last(f["Stress"], step_index))
        E = _tensor(_last(f["Strain"], step_index))
        T = _last(f["Temperature"], step_index).ravel() if "Temperature" in f else None
    if temperature != "field":
        module, func = temperature.rsplit(".", 1)
        T = np.asarray(getattr(_import(module, case_dir), func)(X.T), dtype=float)
        T = np.broadcast_to(T, (len(X),)).copy()
    return {"X": X, "cells": cells, "u": u, "T": T, "S": S, "E": E}


def load_material(path):
    """Card values; E may be a module.func hook, and a card with an eigenstrain
    hook has no alpha (T-dependent card, 'vp1' only)."""
    with open(path) as fh:
        m = yaml.safe_load(fh)
    out = {k: m[k] if isinstance(m[k], str) else float(m[k]) for k in ("E", "nu", "alpha", "T_ref") if k in m}
    for k in ("k", "eigenstrain", "expansion_law"):
        out[k] = m.get(k)
    return out


def _hook(path):
    """Function of a card hook, module.func of the z3st package (materials.steel_1515Ti_Tdep.E_Sobolev)."""
    module, func = path.rsplit(".", 1)
    # import z3st applies its plot style (utils/plotstyle): rcParams, and the titles as no-ops
    # unless Z3ST_PLOT_TITLES is set; the hex_post figures keep their style and titles
    titles = os.environ.get("Z3ST_PLOT_TITLES")
    os.environ["Z3ST_PLOT_TITLES"] = "1"
    try:
        with matplotlib.rc_context():
            return getattr(importlib.import_module(f"z3st.{module}"), func)
    finally:
        if titles is None:
            del os.environ["Z3ST_PLOT_TITLES"]
        else:
            os.environ["Z3ST_PLOT_TITLES"] = titles


# --. geometry and local frame --..
def von_mises(S):
    dev = S - np.trace(S, axis1=1, axis2=2)[:, None, None] / 3 * np.eye(3)
    return np.sqrt(1.5 * np.einsum("cij,cij->c", dev, dev))


def frame(d):
    """Hexagon geometry from the mesh, and the local frame of nodes and cells."""
    X, cells = d["X"], d["cells"]
    r = np.hypot(X[:, 0], X[:, 1])
    corner = X[np.argmax(r)]
    a0 = np.arctan2(corner[1], corner[0]) % (np.pi / 3)

    def local(p):
        sector = np.floor(((np.arctan2(p[:, 1], p[:, 0]) - a0) % (2 * np.pi)) / (np.pi / 3)).astype(int) % 6
        phi_n = a0 + sector * np.pi / 3 + np.pi / 6
        n = np.stack([np.cos(phi_n), np.sin(phi_n), 0 * phi_n], axis=1)
        t = np.stack([-np.sin(phi_n), np.cos(phi_n), 0 * phi_n], axis=1)
        a = (p * n).sum(1)
        return {"sector": sector, "n": n, "t": t, "a": a, "s": (p * t).sum(1), "z": p[:, 2]}

    nodes, cent = local(X), local(X[cells].mean(axis=1))
    a_in, a_out = nodes["a"].min(), nodes["a"].max()
    a_mid = 0.5 * (a_in + a_out)
    L_mid = 2 * a_mid * np.tan(np.pi / 6)
    for f in (nodes, cent):
        f["xi"] = f["a"] - a_mid
        f["s_frac"] = f["s"] / (f["a"] * np.tan(np.pi / 6))
        f["p"] = (f["sector"] + 0.5 * (f["s_frac"] + 1)) * L_mid
    names = []
    for m in range(6):
        deg = int(round(np.degrees(a0 + m * np.pi / 3 + np.pi / 6))) % 360
        names.append(COMPASS.get(deg, f"normal {deg}°"))
    z = X[:, 2]
    return {"nodes": nodes, "cells": cent, "a0": a0, "t": a_out - a_in, "a_mid": a_mid,
            "L_mid": L_mid, "D": 2 * r.max(), "z0": z.min(), "H": np.ptp(z),
            "is_3d": cells.shape[1] == 8, "flats": names}


def perimeter_point(fr, p, xi=0.0):
    """x, y and outward normal of the perimeter coordinate p at wall position xi."""
    m = np.floor(p / fr["L_mid"]).astype(int) % 6
    s_frac = 2 * (p / fr["L_mid"] - np.floor(p / fr["L_mid"])) - 1
    phi_n = fr["a0"] + m * np.pi / 3 + np.pi / 6
    a = fr["a_mid"] + xi
    s = s_frac * a * np.tan(np.pi / 6)
    n = np.stack([np.cos(phi_n), np.sin(phi_n)], axis=-1)
    return a * n[:, 0] - s * n[:, 1], a * n[:, 1] + s * n[:, 0], n


# --. section integrals (2x2x2 Gauss on trilinear hexahedra, VTK order) --..
_REF = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                 [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], dtype=float)
_G = 0.5 + np.array([-1, 1]) * 0.5 / np.sqrt(3)
_GP = np.array([[a, b, c] for a in _G for b in _G for c in _G])


def _shape():
    xi = _GP[:, None, :]
    lin = np.where(_REF[None] == 1, xi, 1 - xi)            # (g, n, 3)
    N = lin.prod(axis=2)
    dN = np.empty((len(_GP), 8, 3))
    sgn = np.where(_REF == 1, 1.0, -1.0)
    for k in range(3):
        others = [j for j in range(3) if j != k]
        dN[:, :, k] = sgn[None, :, k] * lin[:, :, others].prod(axis=2)
    return N, dN


def sections(d, fr):
    """Per z-layer: height, inertia tensor and int eps_zz dx, dy (dx, dy from
    the layer centroid)."""
    X, cells = d["X"], d["cells"]
    N, dN = _shape()
    Xc = X[cells]
    J = np.einsum("gnk,cnd->cgdk", dN, Xc)
    w = np.abs(np.linalg.det(J)) / len(_GP)
    xg = np.einsum("gn,cnd->cgd", N, Xc)

    layers, inv = np.unique(np.round(fr["cells"]["z"], 9), return_inverse=True)
    dz = np.bincount(inv, weights=np.ptp(Xc[:, :, 2], axis=1)) / np.bincount(inv)
    L = lambda v: np.bincount(inv, weights=v) / dz
    A = L(w.sum(1))
    cx, cy = L((w * xg[..., 0]).sum(1)) / A, L((w * xg[..., 1]).sum(1)) / A
    dx = xg[..., 0] - cx[inv][:, None]
    dy = xg[..., 1] - cy[inv][:, None]
    I = np.empty((len(layers), 2, 2))
    I[:, 0, 0], I[:, 1, 1] = L((w * dx * dx).sum(1)), L((w * dy * dy).sum(1))
    I[:, 0, 1] = I[:, 1, 0] = L((w * dx * dy).sum(1))
    ezz = d["E"][:, 2, 2][:, None]
    m_eps = np.stack([L((w * ezz * dx).sum(1)), L((w * ezz * dy).sum(1))], axis=1)
    return {"z": layers, "I": I, "m_eps": m_eps}


def _integrate_twice(z, k):
    """u'' = k along z with u = u' = 0 at z[0]."""
    slope = np.concatenate([[0.0], np.cumsum(0.5 * (k[1:] + k[:-1]) * np.diff(z))])
    return np.concatenate([[0.0], np.cumsum(0.5 * (slope[1:] + slope[:-1]) * np.diff(z))])


def chord(z, u):
    """u measured from the straight line through its two ends."""
    return u - (u[0] + (u[-1] - u[0]) * (z - z[0]) / (z[-1] - z[0]))


def section_mean_u(f, u):
    """Mean displacement of the nodes of each section, at every node (in-plane only)."""
    zu, inv = np.unique(np.round(f["z"], 9), return_inverse=True)
    cnt = np.bincount(inv)
    out = np.zeros_like(u)
    for i in range(2):
        out[:, i] = (np.bincount(inv, weights=u[:, i]) / cnt)[inv]
    return out


def along_z_mean(f, values):
    zu, inv = np.unique(np.round(f["z"], 9), return_inverse=True)
    cnt = np.bincount(inv)
    return zu, np.stack([np.bincount(inv, weights=v) / cnt for v in np.atleast_2d(values.T)], axis=1)


# --. closed forms --..
class Radial:
    """Through-wall gradient: inner flux q, outer T_o, k(T) linear. No bow;
    along s the thick cylinder of the same perimeter and thickness, and the
    slab sigma_ss in the middle half of each flat."""

    bows = False

    def __init__(self, mat, bcs, case_dir, fr):
        th = {bc["type"]: bc for bc in bcs["thermal"]["steel"]}
        self.q = -float(th["Neumann"]["flux"])              # (W/m2) into the wall
        self.T_o = float(th["Dirichlet"]["temperature"])
        module, func = mat["k"].split(".")
        k = getattr(_import(module, case_dir), func)
        self.k0, self.k1 = k(self.T_o), k(self.T_o + 1.0) - k(self.T_o)
        self.mat, self.fr = mat, fr
        self.t = fr["t"]
        self.Rc = 3 * fr["L_mid"] / np.pi                   # same mid-wall perimeter 6 L
        self.a, self.b = self.Rc - self.t / 2, self.Rc + self.t / 2
        self.plane_strain = not fr["is_3d"]

    def _T(self, phi):
        """Kirchhoff inverse: int_{T_o}^{T} k = phi, k linear."""
        if abs(self.k1) < 1e-14:
            return self.T_o + phi / self.k0
        return self.T_o + (-self.k0 + np.sqrt(self.k0**2 + 2 * self.k1 * phi)) / self.k1

    def cylinder(self, xi):
        """Thick cylinder at r = Rc + xi (Timoshenko & Goodier, sec. 150)."""
        E, nu, al, Tr = (self.mat[k] for k in ("E", "nu", "alpha", "T_ref"))
        a, b = self.a, self.b
        rr = np.linspace(a, b, 4001)
        th = self._T(self.q * a * np.log(b / rr)) - Tr
        Ir = np.concatenate([[0.0], np.cumsum(0.5 * (th[1:] * rr[1:] + th[:-1] * rr[:-1]) * np.diff(rr))])
        Ib = Ir[-1]
        r = self.Rc + np.atleast_1d(xi)
        thr, Irr = np.interp(r, rr, th), np.interp(r, rr, Ir)
        c = E * al / (1 - nu)
        s_r = c / r**2 * ((r**2 - a**2) / (b**2 - a**2) * Ib - Irr)
        s_t = c / r**2 * ((r**2 + a**2) / (b**2 - a**2) * Ib + Irr - thr * r**2)
        s_z = c * ((2 * nu if self.plane_strain else 2.0) / (b**2 - a**2) * Ib - thr)
        e_r = (s_r - nu * (s_t + s_z)) / E + al * thr
        e_t = (s_t - nu * (s_r + s_z)) / E + al * thr
        e_z = (s_z - nu * (s_r + s_t)) / E + al * thr
        vm = np.sqrt(0.5 * ((s_r - s_t)**2 + (s_t - s_z)**2 + (s_z - s_r)**2))
        return {"vm": vm, "ss": s_t, "eps_ss": e_t, "eps_zz": e_z, "eps_nn": e_r, "u_n": r * e_t}

    def slab_ss(self, xi):
        """Slab: sigma_ss = alpha E / (1 - nu) (T_mean - T(xi))."""
        E, nu, al = self.mat["E"], self.mat["nu"], self.mat["alpha"]
        x = np.linspace(-self.t / 2, self.t / 2, 2001)
        T = self._T(self.q * (self.t / 2 - x))
        return al * E / (1 - nu) * (T.mean() - np.interp(xi, x, T))

    def references(self, fr, z, layers):
        """Along-s reference curves, per panel key: list of (label, style, p, v)."""
        if fr["is_3d"] and min(z - fr["z0"], fr["z0"] + fr["H"] - z) < 0.2 * fr["H"]:
            return {}, "closed forms valid far from the ends only: FEM only at this height"
        p = np.linspace(0, 6 * fr["L_mid"], 601)
        out = {}
        for side, (xi_c, xi_n) in layers.items():
            cyl_c, cyl_n = self.cylinder(xi_c), self.cylinder(xi_n)
            for key in ("vm", "ss", "eps_ss", "eps_zz", "eps_nn"):
                out.setdefault(key, []).append((f"cylinder, {side}", side, p, np.full_like(p, cyl_c[key][0])))
            out.setdefault("u_n", []).append((f"cylinder, {side}", side, p, np.full_like(p, cyl_n["u_n"][0])))
            window = np.abs(2 * (p / fr["L_mid"] % 1) - 1) < 0.5
            out["ss"].append((f"slab (|s| < L/4), {side}", "slab", p,
                              np.where(window, self.slab_ss(xi_c), np.nan)))
        return out, (f"closed forms: cylinder R = {self.Rc * 1e3:.2f} mm (same perimeter), "
                     f"t = {self.t * 1e3:.2f} mm, same q and T_o; slab = flat far from the vertices")

    def kappa(self, z):
        return np.zeros_like(z)

    def w(self, z):
        return np.zeros_like(z)


class VP1:
    """T linear across the section, from T1 (cold corner) to T2 (hot corner),
    as in IAEA VP1 / VTB hex_duct_linear."""

    bows = True

    def __init__(self, module, mat, case_dir, fr):
        self.m = _import(module, case_dir)
        self.mat, self.fr = mat, fr
        self.e = -np.array([np.cos(self.m.PHI_HOT), np.sin(self.m.PHI_HOT)])   # towards the cold corner
        self.Rc = 3 * fr["L_mid"] / np.pi
        self._zz = np.linspace(fr["z0"], fr["z0"] + fr["H"], 40001)

    def kappa(self, z):
        return self.mat["alpha"] * (self.m.T2(z) - self.m.T1(z)) / self.m.D

    def w(self, z):
        return np.interp(z, self._zz, _integrate_twice(self._zz, self.kappa(self._zz)))

    def u_z(self, z):
        zz = self._zz
        th = 0.5 * (self.m.T1(zz) + self.m.T2(zz)) - self.mat["T_ref"]
        uz = np.concatenate([[0.0], np.cumsum(0.5 * (th[1:] + th[:-1]) * np.diff(zz))]) * self.mat["alpha"]
        return np.interp(z, zz, uz)

    def _linear_at(self, z):
        """T linear in (x, y, z) around z: T1 and T2 constant along z there."""
        return all(abs(f(z + 1e-3) - f(z - 1e-3)) < 1e-9 for f in (self.m.T1, self.m.T2))

    def _free(self, x, y, z):
        """Stress-free solution of a linear T: theta and u at (x, y, z)."""
        T = lambda X, Y: float(np.ravel(self.m.T_vp1(np.array([[X], [Y], [z]])))[0])
        h = 1e-3
        th0 = T(0.0, 0.0) - self.mat["T_ref"]
        g = np.array([(T(h, 0.0) - T(-h, 0.0)) / (2 * h), (T(0.0, h) - T(0.0, -h)) / (2 * h), 0.0])
        P = np.stack([x, y, np.full_like(x, z)], axis=1)
        th = th0 + P @ g
        u = self.mat["alpha"] * (th[:, None] * P - 0.5 * (P**2).sum(1)[:, None] * g)
        return th, u

    def references(self, fr, z, layers):
        if not self._linear_at(z) or z - fr["z0"] < 0.2 * fr["H"]:
            return {}, "no closed form at this height (clamp zone or axial ramp): FEM only"
        al = self.mat["alpha"]
        p = np.linspace(0, 6 * fr["L_mid"], 601)
        x0, y0, _ = perimeter_point(fr, p, 0.0)
        ang = np.arctan2(y0, x0)                              # cylinder at the same polar angle
        out = {}
        for side, (xi_c, xi_n) in layers.items():
            x, y, _ = perimeter_point(fr, p, xi_c)
            th, _ = self._free(x, y, z)
            for key, v in (("vm", 0 * th), ("ss", 0 * th), ("eps_ss", al * th),
                           ("eps_zz", al * th), ("eps_nn", al * th)):
                out.setdefault(key, []).append((f"hexagon closed form, {side}", "hex", p, v))
            x, y, n = perimeter_point(fr, p, xi_n)
            _, u = self._free(x, y, z)
            u0 = self._free(np.zeros(1), np.zeros(1), z)[1][0, :2]          # section centre
            out.setdefault("u_n", []).append((f"hexagon closed form, {side}", "hex", p,
                                              ((u[:, :2] - u0) * n).sum(1)))
            for xi, keys in ((xi_c, ("vm", "ss", "eps_ss", "eps_zz", "eps_nn")), (xi_n, ("u_n",))):
                r = self.Rc + xi
                th, u = self._free(r * np.cos(ang), r * np.sin(ang), z)
                vals = {"vm": 0 * th, "ss": 0 * th, "eps_ss": al * th, "eps_zz": al * th, "eps_nn": al * th,
                        "u_n": (u[:, 0] - u0[0]) * np.cos(ang) + (u[:, 1] - u0[1]) * np.sin(ang)}
                for key in keys:
                    out[key].append((f"cylinder, {side}", side, p, vals[key]))
        return out, (f"closed forms (T linear: stress-free far from the clamp, rigid motion of the base zone "
                     f"not included); cylinder R = {self.Rc * 1e3:.2f} mm (same perimeter) at the same polar angle")


class VP1Tdep(VP1):
    """VP1 with a T-dependent card (E(T) hook, eigenstrain hook in place of
    alpha), as VP1/1515Ti_Tdependent/non-regression.py: on each section
    eps_th(T) of the mid-wall line is fitted by plane sections, eps ~ a + b.(x, y),
    weighted by E(T) ds; kappa = -b.e, u_z = int a dz. For constant alpha and E
    it is the VP1 closed form. Thermal part of the hook only (no swelling, no
    imposed creep in VP1)."""

    def __init__(self, module, mat, case_dir, fr):
        super().__init__(module, mat, case_dir, fr)
        E = _hook(mat["E"])
        eps_th = _hook(mat["eigenstrain"].rsplit(".", 1)[0] + ".thermal_strain")
        law = mat["expansion_law"] or "Sobolev"
        p = (np.arange(6 * 400) + 0.5) * fr["L_mid"] / 400       # mid-wall line, 400 points a flat
        x, y, _ = perimeter_point(fr, p, 0.0)
        A = np.stack([np.ones_like(x), x, y], axis=1)
        self._zk = np.linspace(fr["z0"], fr["z0"] + fr["H"], 4001)
        c = np.empty((len(self._zk), 3))
        for i, z in enumerate(self._zk):
            T = self.m.T_vp1(np.array([x, y, np.full_like(x, z)]))
            w = E(T)
            c[i] = np.linalg.solve(A.T @ (w[:, None] * A), A.T @ (w * eps_th(T, mat["T_ref"], law)))
        self._a, self._k = c[:, 0], -c[:, 1:] @ self.e

    def kappa(self, z):
        return np.interp(z, self._zk, self._k)

    def u_z(self, z):
        uz = np.concatenate([[0.0], np.cumsum(0.5 * (self._a[1:] + self._a[:-1]) * np.diff(self._zk))])
        return np.interp(z, self._zk, uz)

    def references(self, fr, z, layers):
        return {}, "T-dependent card: eps_th(T) not linear in x, no stress-free solution along s: FEM only"


def closed_form(cfg, mat, case_dir, fr):
    kind = str(cfg.get("closed_form", "none"))
    if kind == "radial":
        with open(os.path.join(case_dir, cfg["boundary_conditions"])) as fh:
            return Radial(mat, yaml.safe_load(fh), case_dir, fr)
    if kind.startswith("vp1:"):
        mod_dir, module = os.path.split(os.path.normpath(os.path.join(case_dir, kind.split(":", 1)[1])))
        return (VP1 if "alpha" in mat else VP1Tdep)(module, mat, mod_dir, fr)
    return None


# --. analysis --..
def analyse(d, mat, support="free", closed=None):
    """FEM quantities of one dataset, projected on the closed-form bowing
    direction when there is one."""
    fr = frame(d)
    c = fr["cells"]
    S, E = d["S"], d["E"]
    res = {"frame": fr, "support": support, "data": d, "closed": closed,
           "vm": von_mises(S),
           "ss": np.einsum("ci,cij,cj->c", c["t"], S, c["t"]),
           "eps_ss": np.einsum("ci,cij,cj->c", c["t"], E, c["t"]),
           "eps_nn": np.einsum("ci,cij,cj->c", c["n"], E, c["n"]),
           "eps_zz": E[:, 2, 2],
           "u_n": ((d["u"] - section_mean_u(fr["nodes"], d["u"])) * fr["nodes"]["n"]).sum(1)}
    if not fr["is_3d"]:
        return res
    sec = sections(d, fr)
    kap = -np.einsum("lij,lj->li", np.linalg.inv(sec["I"]), sec["m_eps"])
    if closed is not None and closed.bows:
        e = closed.e
    else:
        k_max = kap[np.argmax(np.linalg.norm(kap, axis=1))]
        e = k_max / np.linalg.norm(k_max) if np.linalg.norm(k_max) > 1e-9 else np.array([1.0, 0.0])
    zn, ubar = along_z_mean(fr["nodes"], d["u"])
    w = ubar[:, :2] @ e
    res.update({"e": e, "e_perp": np.array([-e[1], e[0]]), "z_sec": sec["z"], "kappa": kap @ e,
                "z_nodes": zn, "u_bar": ubar, "w": w if support == "clamped" else chord(zn, w)})
    return res


def max_stresses(res):
    fr, c = res["frame"], res["frame"]["cells"]
    C = res["data"]["X"][res["data"]["cells"]].mean(axis=1)
    szz = res["data"]["S"][:, 2, 2]
    out = {}
    for name, v, k in (("von Mises", res["vm"], np.argmax(res["vm"])),
                       ("sigma_ss max", res["ss"], np.argmax(res["ss"])),
                       ("sigma_ss min", res["ss"], np.argmin(res["ss"])),
                       ("|sigma_zz| max", szz, np.argmax(np.abs(szz)))):
        p = c["p"][k] % fr["L_mid"]
        out[name] = {"MPa": float(v[k] * 1e-6), "x_m": float(C[k, 0]), "y_m": float(C[k, 1]),
                     "z_m": float(C[k, 2]), "flat": fr["flats"][c["sector"][k]],
                     "s_mm": float(c["s"][k] * 1e3), "xi_mm": float(c["xi"][k] * 1e3),
                     "to_vertex_mm": float(min(p, fr["L_mid"] - p) * 1e3)}
    return out


def stations(res, cfg_z=None):
    fr = res["frame"]
    if cfg_z:
        return np.array(cfg_z, dtype=float)
    return fr["z0"] + fr["H"] * np.array([0.25, 0.5, 0.75, 1.0])


def beam_table(runs, zs):
    """Closed form and every run at the stations zs; a run with a closed form
    other than that of the first run also gets its own closed-form values."""
    cf = runs[0][1]["closed"]
    rows = []
    for z in zs:
        row = {"z_m": float(z)}
        if cf is not None:
            row["w_closed_mm"] = float(cf.w(np.array([z]))[0] * 1e3)
            row["kappa_closed"] = float(cf.kappa(np.array([z]))[0])
        for label, r in runs:
            row[f"w_{label}_mm"] = float(np.interp(z, r["z_nodes"], r["w"]) * 1e3)
            row[f"kappa_{label}"] = float(np.interp(z, r["z_sec"], r["kappa"]))
            if r["closed"] is not None and r["closed"] is not cf:
                row[f"w_closed_{label}_mm"] = float(r["closed"].w(np.array([z]))[0] * 1e3)
                row[f"kappa_closed_{label}"] = float(r["closed"].kappa(np.array([z]))[0])
        rows.append(row)
    return rows


def _closed_forms(runs):
    """Distinct closed forms of the runs: (label of the first run with it, closed form, colour)."""
    out = []
    for (label, r), color in zip(runs, RUN_COLORS):
        if r["closed"] is not None and all(r["closed"] is not c for _, c, _ in out):
            out.append((label, r["closed"], color))
    return out


# --. figures --..
def _finish(fig, title, subtitle, top=0.90):
    fig.suptitle(title, fontsize=15, fontweight="bold", y=0.995)
    fig.text(0.5, 0.96, subtitle, ha="center", va="top", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, top))
    return fig


def _table(ax, header, rows, fill=False):
    """Table on ``ax``; ``fill`` stretches it over the whole axes."""
    ax.axis("off")
    kw = {"bbox": [0, 0, 1, 1]} if fill else {"loc": "center"}
    tb = ax.table(cellText=rows, colLabels=header, cellLoc="center", **kw)
    tb.auto_set_font_size(False)
    tb.set_fontsize(9)
    if not fill:
        tb.scale(1.0, 1.6)
    for (r, _), cell in tb.get_celld().items():
        if r == 0:
            cell.set_facecolor("#e8e8e8")
            cell.set_text_props(fontweight="bold")
    return tb


def _zoom_bulk(ax, data, fr, refs=(), ratio=5.0):
    """Full range on ``ax``; when the ``data`` within one D of the bottom and
    of the top (Saint-Venant end zones) stretch the x range to more than
    ``ratio`` times that of the rest, an inset zoomed on the rest. ``data``,
    ``refs``: (x, z_mm, style) with style the plot kwargs; ``refs`` (reference
    lines) are drawn in the inset only."""
    z_lo, z_hi = (fr["z0"] + fr["D"]) * 1e3, (fr["z0"] + fr["H"] - fr["D"]) * 1e3
    allx = np.concatenate([x for x, _, _ in data])
    bulk = np.concatenate([x[(z > z_lo) & (z < z_hi)] for x, z, _ in data])
    if not bulk.size or np.ptp(allx) <= ratio * max(np.ptp(bulk), 1e-12 * np.abs(allx).max()):
        return
    lo, hi = bulk.min(), bulk.max()
    pad = 0.1 * max(hi - lo, 1e-3 * np.abs(bulk).max(), 1e-30)
    x0, x1 = ax.get_xlim()
    left = 0.5 * (lo + hi) > 0.5 * (x0 + x1)     # inset on the side away from the bulk
    ins = ax.inset_axes([0.06 if left else 0.5, 0.25, 0.44, 0.45])
    for x, z, style in list(data) + list(refs):
        ins.plot(x, z, **style)
    ins.set_xlim(lo - pad, hi + pad)
    ins.set_ylim(ax.get_ylim())
    ins.set_title(f"zoom: farther than D = {fr['D'] * 1e3:.0f} mm from the ends", fontsize=8)
    ins.tick_params(labelsize=7)
    ins.grid(alpha=0.4)


def fig_beam(runs, zs, title, kind):
    """VTB layout: quantity on x, axial position (mm) on y; closed form as
    circles on the same axes; table at the stations."""
    r0 = runs[0][1]
    cf = r0["closed"]
    is_w = kind == "w"
    zkey = "z_nodes" if is_w else "z_sec"
    scale, unit = (1e3, "mm") if is_w else (1.0, "1/m")
    name = "Bowing displacement" if is_w else "Curvature κ"
    rel = cf is not None and cf.bows
    if rel:   # curves, their difference from the closed form, table below
        fig = plt.figure(figsize=(14, 11))
        gs = fig.add_gridspec(2, 2, height_ratios=(3, 1.1), hspace=0.3)
        ax, tax = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, :])
        dax = fig.add_subplot(gs[0, 1], sharey=ax)
    else:     # curve (or note) next to the table
        fig = plt.figure(figsize=(14, 6.5))
        gs = fig.add_gridspec(1, 2, width_ratios=(1, 1.1), wspace=0.15)
        ax, tax, dax = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), None

    fr = r0["frame"]
    fem_max = max(np.abs(r[kind]).max() for _, r in runs)
    cfs = _closed_forms(runs)
    own = lambda r: r["closed"] if r["closed"] is not None else cf     # closed form a run is compared with
    if rel:   # one difference for the panel and the table, each run against its own closed form
        if is_w:
            dunit, diff = "µm", (lambda v, ref, c: (v - ref) * 1e6)
        else:     # in % of the largest closed-form κ: defined also where κ_closed = 0
            zg = np.linspace(fr["z0"], fr["z0"] + fr["H"], 401)
            dunit, diff = "% of max closed-form κ", (lambda v, ref, c: (v - ref) / np.abs(c.kappa(zg)).max() * 100)
    if cf is None or cf.bows:
        lines = [(r[kind] * scale, r[zkey] * 1e3, dict(color=color, lw=1.6, label=f"FEM {label}"))
                 for (label, r), color in zip(runs, RUN_COLORS)]
        points = []
        zc = np.linspace(fr["z0"], fr["z0"] + fr["H"], 26)
        for lab, c, color in cfs:    # red circles; with several closed forms, the colour of their run
            v = c.w(zc) if is_w else c.kappa(zc)
            one = len(cfs) == 1
            points.append((v * scale, zc * 1e3, dict(ls="none", marker="o", mfc="none", ms=7,
                                                     mec="red" if one else color,
                                                     label="closed form (equation)" if one else f"closed form, {lab}")))
        for x, z, style in lines + points:
            ax.plot(x, z, **style)
        ax.set(xlabel=f"{name} ({unit})", ylabel="Axial location (mm)")
        ax.grid(alpha=0.4)
        ax.legend(loc="upper left", fontsize=9)
        _zoom_bulk(ax, lines + points, fr)
    if rel:   # difference panel
        lines = [(diff(r[kind], own(r).w(r[zkey]) if is_w else own(r).kappa(r[zkey]), own(r)), r[zkey] * 1e3,
                  dict(color=color, lw=1.6, label=f"FEM {label}")) for (label, r), color in zip(runs, RUN_COLORS)]
        for x, z, style in lines:
            dax.plot(x, z, **style)
        dax.axvline(0, color="red", lw=1.0)
        dax.set(xlabel=f"FEM − closed form ({dunit})", title="Difference from the closed form")
        dax.grid(alpha=0.4)
        dax.legend(loc="upper left", fontsize=9)
        plt.setp(dax.get_yticklabels(), visible=False)
        _zoom_bulk(dax, lines, fr, refs=[(np.zeros(2), np.array(ax.get_ylim()), dict(color="red", lw=1.0))])
    if not (cf is None or cf.bows):
        ax.axis("off")
        ax.text(0.5, 0.5, f"Closed form: {name.lower()} = 0\n(temperature symmetric under\nthe hexagon rotations)\n\n"
                f"FEM: max |{'w' if is_w else 'κ'}| = {fem_max * scale:.1e} {unit}\n(round-off, nothing to plot)",
                ha="center", va="center", fontsize=12, transform=ax.transAxes,
                bbox=dict(boxstyle="round", fc="#f4f4f4", ec="0.7"))

    rows = beam_table(runs, zs)
    key = (lambda lab: f"w_{lab}_mm") if is_w else (lambda lab: f"kappa_{lab}")
    ck = "w_closed_mm" if is_w else "kappa_closed"
    ckey = lambda lab, c: ck if c is cf else (f"w_closed_{lab}_mm" if is_w else f"kappa_closed_{lab}")
    fmt = (lambda v: f"{v + 0.0:.3f}".replace("-0.000", "0.000")) if is_w else (lambda v: f"{v:.4e}")
    header = (["z (mm)"] + (["closed form"] if len(cfs) == 1 else [f"closed form\n{lab}" for lab, _, _ in cfs])
              + [f"FEM\n{lab}" for lab, _ in runs] + ([f"{lab}\n− closed form" for lab, _ in runs] if rel else []))
    si = 1e-3 if is_w else 1.0     # table w in mm
    body = []
    for row in rows:
        line = [f"{row['z_m'] * 1e3:.0f}"] + [fmt(row[ckey(lab, c)]) for lab, c, _ in cfs]
        line += [fmt(row[key(lab)]) for lab, _ in runs]
        if rel:
            dv = [f"{diff(row[key(lab)] * si, row[ckey(lab, own(r))] * si, own(r)):+.2f}" for lab, r in runs]
            line += [s.replace("-0.00", "0.00").replace("+0.00", "0.00") for s in dv]
        body.append(line)
    if rel:
        tb = _table(tax, header, body, fill=True)
        for (r, _), cell in tb.get_celld().items():
            if r == 0:
                cell.set_height(cell.get_height() * 1.6)
    else:
        _table(tax, header, body)
    ttl = f"{name} ({unit}) at the stations"
    if rel:
        ttl += f"; differences from the closed form in {dunit}"
    tax.set_title(ttl, fontsize=11)
    sub = f"{name} along the axis"
    if cf is None:
        sub += " — no closed form for this case: FEM only"
    elif not cf.bows:
        sub += " — closed form zero by symmetry: values in the table"
    else:
        sub += " — FEM (lines) vs closed form (circles)"
    if not is_w:
        sub += "\nFEM curvature: plane-section fit of ε_zz over each section"
    elif r0["support"] == "clamped":
        sub += "\nclamped base, u = u' = 0; FEM: mean of the nodes of each section"
    else:
        sub += "\nfree ends: measured from the chord between the two ends; FEM: mean of the nodes of each section"
    fig.suptitle(title, fontsize=15, fontweight="bold", y=0.995)
    fig.text(0.5, 0.965, sub, ha="center", va="top", fontsize=10)
    fig.subplots_adjust(top=0.89 if rel else 0.83, bottom=0.04 if rel else 0.09, left=0.07, right=0.97)
    return fig


def fig_displacement(runs, title):
    r0 = runs[0][1]
    cf = r0["closed"]
    e, ep = r0["e"], r0["e_perp"]
    fig, ax = plt.subplots(1, 3, figsize=(15, 7), sharey=True)
    comps = (("u along the bowing direction", lambda r: r["u_bar"][:, :2] @ e),
             ("u perpendicular to it", lambda r: r["u_bar"][:, :2] @ ep),
             ("u_z", lambda r: r["u_bar"][:, 2]))
    for i, (name, get) in enumerate(comps):
        for (label, r), color in zip(runs, RUN_COLORS):
            ax[i].plot(get(r) * 1e3, r["z_nodes"] * 1e3, "-", color=color, lw=1.6, label=f"FEM {label}")
        ax[i].set(xlabel=f"section mean {name} (mm)")
        ax[i].grid(alpha=0.4)
    zc = np.linspace(r0["frame"]["z0"], r0["frame"]["z0"] + r0["frame"]["H"], 26)
    ref = dict(ls="none", marker="o", mfc="none", mec="red", label="closed form")
    if cf is not None and cf.bows and r0["support"] == "clamped":
        cfs = _closed_forms(runs)
        for lab, c, color in cfs:    # with several closed forms, the colour of their run
            style = ref if len(cfs) == 1 else dict(ref, mec=color, label=f"closed form, {lab}")
            ax[0].plot(c.w(zc) * 1e3, zc * 1e3, **style)
            ax[1].plot(0 * zc, zc * 1e3, **style)
            ax[2].plot(c.u_z(zc) * 1e3, zc * 1e3, **style)
    elif cf is not None and not cf.bows:
        for a in ax[:2]:
            a.plot(0 * zc, zc * 1e3, **ref)
    ax[0].set_ylabel("Axial location (mm)")
    for a in ax:
        a.legend(fontsize=8)
    sub = "Displacement along the axis (mean of the nodes of each section)"
    if cf is None:
        sub += " — no closed form: FEM only"
    if r0["support"] != "clamped":
        sub += "\nrigid-body modes removed by the solver: u_x, u_y defined up to a rigid motion"
    return _finish(fig, title, sub, top=0.91)


def _max_per_z(z, v):
    zu, inv = np.unique(np.round(z, 9), return_inverse=True)
    out = np.full(len(zu), -np.inf)
    np.maximum.at(out, inv, v)
    return zu, out


def fig_stress_z(runs, title):
    fig, ax = plt.subplots(1, 2, figsize=(13, 7), sharey=True)
    for (label, r), color in zip(runs, RUN_COLORS):
        zc = r["frame"]["cells"]["z"]
        z, v = _max_per_z(zc, r["vm"] * 1e-6)
        ax[0].plot(v, z * 1e3, "-", color=color, label=f"FEM {label}")
        z, v = _max_per_z(zc, np.abs(r["ss"]) * 1e-6)
        ax[1].plot(v, z * 1e3, "-", color=color, label=f"FEM {label}")
        m = max_stresses(r)["von Mises"]
        ax[0].plot(m["MPa"], m["z_m"] * 1e3, "*", color=color, ms=14, mec="k",
                   label=f"max {label}: {m['MPa']:.1f} MPa at z = {m['z_m'] * 1e3:.1f} mm")
    ax[0].set(xlabel="max σ_vm in the section (MPa)", ylabel="Axial location (mm)")
    ax[1].set(xlabel="max |σ_ss| in the section (MPa)")
    for a in ax:
        a.set_xscale("symlog", linthresh=1.0)
        a.grid(alpha=0.4)
        a.legend(fontsize=8)
    return _finish(fig, title, "Largest stresses of each section along the axis (linear below 1 MPa, log above)",
                   top=0.93)


def _nearest(v, target):
    d = np.abs(v - target)
    return d < d.min() + 1e-9


def _ring(f, field, z, xi):
    m = _nearest(f["z"], z) & _nearest(f["xi"], xi)
    pu, inv = np.unique(np.round(f["p"][m], 9), return_inverse=True)
    return pu, np.bincount(inv, weights=field[m]) / np.bincount(inv)


def _layers(fr):
    """Wall layers drawn along s: {side: (cell xi, node xi)}, and whether
    there is a single cell through the wall."""
    xc = np.unique(np.round(fr["cells"]["xi"], 9))
    xn = np.unique(np.round(fr["nodes"]["xi"], 9))
    return {"inner": (xc[0], xn[0]), "outer": (xc[-1], xn[-1])}, len(xc) == 1


PANELS = (("vm", 1e-6, "σ_vm (MPa)"), ("ss", 1e-6, "hoop σ_ss (MPa)"),
          ("eps_ss", 1.0, "ε_ss (-)"), ("eps_zz", 1.0, "ε_zz (-)"),
          ("eps_nn", 1.0, "ε_nn (-)"), ("u_n", 1e3, "u_n relative to the section centre (mm)"))
REF_STYLE = {"inner": dict(color="0.25", ls=":", lw=2.2),
             "outer": dict(color="0.6", ls=":", lw=2.2),
             "hex": dict(color="red", ls="none", marker="o", ms=3.5, mfc="none", markevery=10),
             "slab": dict(color="k", ls="-", lw=4.0, alpha=0.3)}


PANEL = {key: (scale, ylabel) for key, scale, ylabel in PANELS}
# along-s figures at the post.yaml heights (along_s_z): one png per pair of panels
ALONG_S_GROUPS = (("stress", ("ss", "vm")), ("strain_sz", ("eps_ss", "eps_zz")), ("normal", ("eps_nn", "u_n")))
RECAP_LS, RECAP_NAMES = ("-", "--", "-."), ("solid", "dashed", "dash-dot")   # one per run in the recap
LEVEL_COLORS = ("tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown")


def _references(r0, z):
    fr0, cf = r0["frame"], r0["closed"]
    if cf is None:
        return {}, "no closed form for this case: FEM only"
    return cf.references(fr0, z, _layers(fr0)[0])


def _sides(key, single):
    """Wall layers drawn: with one cell through the wall the cell fields have
    a single (mid-wall) layer."""
    return ("inner", "outer") if (key == "u_n" or not single) else ("inner",)


def _draw_s(ax, runs, z, key, refs):
    """Closed forms and FEM of every run along s at z, for one panel key."""
    scale = PANEL[key][0]
    sides = _sides(key, _layers(runs[0][1]["frame"])[1])
    for label, side, p, v in refs.get(key, []):
        if side in ("inner", "outer") and side not in sides:
            continue
        if label.endswith(", outer") and "outer" not in sides:
            continue
        lab = label if len(sides) > 1 else label.rsplit(",", 1)[0] + ", mid-wall"
        ax.plot(p * 1e3, v * scale, label=lab, **REF_STYLE[side])
    for (label, r), color in zip(runs, RUN_COLORS):
        f = r["frame"]["nodes"] if key == "u_n" else r["frame"]["cells"]
        lay, _ = _layers(r["frame"])
        for side in sides:
            xi_c, xi_n = lay[side]
            p, v = _ring(f, r[key], z, xi_n if key == "u_n" else xi_c)
            lc = color if len(runs) > 1 else ("tab:blue" if side == "inner" else "tab:red")
            tag = side if len(sides) > 1 else "mid-wall"
            ax.plot(p * 1e3, v * scale, LAYER_STYLES[side], color=lc, lw=1.5, label=f"FEM {label}, {tag}")


def _s_axes(ax, fr0, key, flats=False, xlabel=False):
    L = fr0["L_mid"] * 1e3
    for m in range(7):
        ax.axvline(m * L, color="k", lw=1.0, alpha=0.55)
    ax.set_ylabel(PANEL[key][1])
    ax.ticklabel_format(axis="y", useOffset=False)
    ax.grid(alpha=0.35)
    ax.legend(fontsize=7, loc="best", ncol=2)
    if flats:
        for m in range(6):
            ax.text((m + 0.5) * L, 1.01, fr0["flats"][m], ha="center", va="bottom", fontsize=9,
                    transform=ax.get_xaxis_transform())
    if xlabel:
        ax.set_xlabel("s along the perimeter at mid-wall (mm); vertical lines = vertices V1..V6")
        ax.set_xticks(np.arange(7) * L)
        ax.set_xticklabels([f"V{m % 6 + 1}\n{m * L:.0f}" for m in range(7)])


def _layer_note(single):
    return ("one cell through the wall: stresses and strains at mid-wall, u_n on the inner/outer faces"
            if single else "cell layers next to the inner (solid) and outer (dashed) faces")


def fig_along_s(runs, z, title):
    r0 = runs[0][1]
    fr0 = r0["frame"]
    refs, ref_note = _references(r0, z)
    fig, axes = plt.subplots(3, 2, figsize=(16, 12.5), sharex=True)
    for i, (ax, (key, _, _)) in enumerate(zip(axes.flat, PANELS)):
        _draw_s(ax, runs, z, key, refs)
        _s_axes(ax, fr0, key, flats=i < 2, xlabel=i >= 4)
    where = f"z = {z * 1e3:.1f} mm" if fr0["is_3d"] else "plane-strain section"
    return _finish(fig, title, f"Along s, {where} — {_layer_note(_layers(fr0)[1])}\n{ref_note}", top=0.92)


def _z_cells(fr, z):
    """Cell layers used for the height z (two when z falls on a cell face)."""
    zc = fr["cells"]["z"]
    return np.unique(np.round(zc[_nearest(zc, z)], 9))


def fig_along_s_levels(runs, levels, keys, title):
    """Two panels (``keys``) along s at every height of ``levels`` {name: z},
    one row each, then a recap row with the heights overlaid."""
    r0 = runs[0][1]
    fr0 = r0["frame"]
    single = _layers(fr0)[1]
    n = len(levels)
    fig, axes = plt.subplots(n + 1, 2, figsize=(16, 3.3 * (n + 1) + 1.2), sharex=True)
    for i, (name, z) in enumerate(levels.items()):
        refs, note = _references(r0, z)
        cells = "/".join(f"{c * 1e3:.0f}" for c in _z_cells(fr0, z))
        ref = "closed form shown" if refs else "FEM only: no closed form at this height"
        for j, key in enumerate(keys):
            ax = axes[i, j]
            _draw_s(ax, runs, z, key, refs)
            _s_axes(ax, fr0, key, flats=i == 0)
            ax.set_title(f"{name}: z = {z * 1e3:.0f} mm (cell layer {cells} mm) — {ref}", fontsize=9,
                         loc="left", pad=16 if i == 0 else 4)
    # recap: color = height; solid first run, dashed the others, dotted the closed form
    for j, key in enumerate(keys):
        ax = axes[n, j]
        scale = PANEL[key][0]
        side = "outer" if key == "u_n" else "inner"
        for (name, z), color in zip(levels.items(), LEVEL_COLORS):
            for k, (label, r) in enumerate(runs):
                f = r["frame"]["nodes"] if key == "u_n" else r["frame"]["cells"]
                xi_c, xi_n = _layers(r["frame"])[0][side]
                p, v = _ring(f, r[key], z, xi_n if key == "u_n" else xi_c)
                ax.plot(p * 1e3, v * scale, RECAP_LS[k % len(RECAP_LS)], color=color, lw=1.5,
                        label=name if k == 0 else None)
            refs, _ = _references(r0, z)
            cand = [c for c in refs.get(key, []) if c[1] == "hex" and (side in c[0] or single)]
            for lab, _, p, v in cand[:1]:
                ax.plot(p * 1e3, v * scale, ":", color=color, lw=2.2)
        _s_axes(ax, fr0, key, xlabel=True)
        wall = f"{side} face" if key == "u_n" else ("mid-wall" if single else "inner layer")
        others = ", ".join(f"{RECAP_NAMES[k % len(RECAP_LS)]} {lab}" for k, (lab, _) in enumerate(runs) if k)
        ax.set_title(f"recap: all heights overlaid ({wall}; solid FEM {runs[0][0]}"
                     + (f", {others}" if others else "") + ", dotted closed form where it exists)",
                     fontsize=9, loc="left")
    sub = (f"Along s at {n} heights (rows) and recap (bottom) — {_layer_note(single)}\n"
           + _references(r0, max(levels.values()))[1])
    return _finish(fig, title, sub, top=1 - 1.3 / (3.3 * (n + 1) + 1.2))


def fig_max_table(runs, title):
    rows = []
    for label, r in runs:
        for name, v in max_stresses(r).items():
            rows.append([label, name, f"{v['MPa']:.2f}", v["flat"], f"{v['s_mm']:+.2f}",
                         f"{v['to_vertex_mm']:.2f}", f"{v['xi_mm']:+.2f}", f"{v['z_m'] * 1e3:.1f}",
                         f"({v['x_m'] * 1e3:+.1f}, {v['y_m'] * 1e3:+.1f})"])
    fig, ax = plt.subplots(figsize=(15, 1.3 + 0.36 * (len(rows) + 1)))
    _table(ax, ["code", "quantity", "value (MPa)", "flat", "s (mm)", "to vertex (mm)", "ξ (mm)",
                "z (mm)", "(x, y) (mm)"], rows, fill=True)
    return _finish(fig, title, "Maximum stresses: magnitude and position (cell values; s from mid-flat, "
                   "ξ from mid-wall, outward)", top=0.85)


def _crop(img, pad=10):
    """Trim the white border of a render."""
    ink = np.where((img[..., :3] < 250).any(axis=2))
    r0, r1 = max(ink[0].min() - pad, 0), min(ink[0].max() + pad, img.shape[0])
    c0, c1 = max(ink[1].min() - pad, 0), min(ink[1].max() + pad, img.shape[1])
    return img[r0:r1, c0:c1]


def fig_vtb_3d(res, title, path):
    """VTB figure: temperature (left) and the deformed duct, magnified,
    coloured by the displacement along the bowing direction (right)."""
    import pyvista as pv
    d, fr = res["data"], res["frame"]
    grid = pv.UnstructuredGrid({pv.CellType.HEXAHEDRON: d["cells"]}, d["X"])
    T_c = d["T"] - 273.15
    grid.point_data["T"] = T_c
    grid.point_data["u"] = d["u"]
    e = res["e"]
    bows = np.abs(res["w"]).max() > 1e-6
    if bows:
        val, label = d["u"][:, :2] @ e * 1e3, "displacement along the bowing direction (mm)"
        amp = 0.05 * fr["H"] / np.abs(res["w"]).max()
    else:
        val, label = np.linalg.norm(d["u"], axis=1) * 1e3, "|u| (mm)"
        amp = 0.05 * fr["D"] / max(val.max() * 1e-3, 1e-30)
    amp = max(1.0, round(amp))
    grid.point_data["v"] = val
    warped = grid.warp_by_vector("u", factor=amp)
    view = np.array([-e[1], e[0], 0.0])                    # looking across the bowing plane
    c = np.array(grid.center)
    panels = []
    for mesh, sc, v, lab, sub in ((grid, "T", T_c, "Temperature (°C)", "Temperature"),
                                  (warped, "v", val, label, f"Deformed shape (x{amp:.0f})")):
        p = pv.Plotter(off_screen=True, window_size=(600, 1600))
        p.set_background("white")
        p.add_mesh(mesh.extract_surface(algorithm="dataset_surface"), scalars=sc, cmap="coolwarm",
                   clim=(v.min(), v.max()), show_scalar_bar=False)
        p.camera_position = [tuple(c - view * 10 * fr["H"]), tuple(c), (0, 0, 1)]   # e to the right
        p.enable_parallel_projection()
        p.reset_camera()
        panels.append((_crop(p.screenshot(return_img=True)), v, lab, sub))
        p.close()
    fig, ax = plt.subplots(1, 2, figsize=(9, 10))
    for a, (img, v, lab, sub) in zip(ax, panels):
        im = a.imshow(img)
        a.axis("off")
        a.set_title(sub, fontsize=11)
        sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=plt.Normalize(v.min(), v.max()))
        fig.colorbar(sm, ax=a, fraction=0.08, shrink=0.55, location="left", pad=0.02, label=lab)
    sub = ("Temperature mapping (left) and duct bowing (right, magnified), as VTB hex_duct_linear"
           if bows else "Temperature mapping (left) and deformed shape (right, magnified): no bowing")
    sub += f"\nview across the bowing plane: the bowing direction points to the right"
    _finish(fig, title, sub, top=0.93).savefig(path, dpi=130)
    plt.close(fig)
    print(f"[INFO] {path}")


def s_heights(res):
    """Heights of the along-s figures: the max von Mises, and mid-height."""
    fr = res["frame"]
    if not fr["is_3d"]:
        return {"section": 0.0}
    zc = fr["cells"]["z"]
    z_vm = zc[np.argmax(res["vm"])]
    z_mid = zc[_nearest(zc, fr["z0"] + fr["H"] / 2)][0]
    out = {"zmaxvm": z_vm}
    if abs(z_mid - z_vm) > 1e-9:
        out["zmid"] = z_mid
    return out


def save_all(runs, out, prefix, title, zs=None, vtb_3d=True, levels=None):
    """Every figure of ``runs`` as ``out/<prefix><figure>.png``, titled ``title``.
    With ``levels`` {name: z} the along-s figures are along_s_<group>.png at
    those heights instead of along_s_<tag>.png at s_heights."""
    os.makedirs(out, exist_ok=True)
    r0 = runs[0][1]
    figs = []
    if r0["frame"]["is_3d"]:
        zs = stations(r0, zs)
        if vtb_3d:
            fig_vtb_3d(r0, title, os.path.join(out, f"{prefix}vtb_temperature_bowing.png"))
        figs += [("vtb_bowing", fig_beam(runs, zs, title, "w")),
                 ("curvature", fig_beam(runs, zs, title, "kappa")),
                 ("displacement_z", fig_displacement(runs, title)),
                 ("stress_z", fig_stress_z(runs, title))]
    if levels:
        for group, keys in ALONG_S_GROUPS:
            figs.append((f"along_s_{group}", fig_along_s_levels(runs, levels, keys, title)))
    else:
        for tag, z in s_heights(r0).items():
            figs.append((f"along_s_{tag}", fig_along_s(runs, z, title)))
    figs.append(("max_stress", fig_max_table(runs, title)))
    for name, fig in figs:
        path = os.path.join(out, f"{prefix}{name}.png")
        fig.savefig(path, dpi=130)
        plt.close(fig)
        print(f"[INFO] {path}")
    summary = {"max_stresses": {lab: max_stresses(r) for lab, r in runs}}
    if r0["frame"]["is_3d"]:
        summary["stations"] = beam_table(runs, zs)
    with open(os.path.join(out, f"{prefix}summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    return summary


def run_case(case_dir):
    case_dir = os.path.abspath(case_dir)
    with open(os.path.join(case_dir, "post.yaml")) as fh:
        cfg = yaml.safe_load(fh)
    path = lambda p: os.path.normpath(os.path.join(case_dir, p))
    xdmf = path(cfg["fields"])
    case = os.path.basename(case_dir)
    if not os.path.isfile(os.path.splitext(xdmf)[0] + ".h5"):
        print(f"[SKIP] {case}: no {os.path.relpath(xdmf, case_dir)} (run the case first)")
        return None
    d = load_z3st(xdmf, temperature=cfg.get("temperature", "field"), case_dir=case_dir)
    if d["T"] is None:
        sys.exit(f"{case}: no Temperature in {xdmf}; set 'temperature: module.func' in post.yaml")
    mat = load_material(path(cfg["material"]))
    cf = closed_form(cfg, mat, case_dir, frame(d))
    res = analyse(d, mat, cfg.get("support", "free"), cf)
    save_all([("z3st", res)], RESULTS, f"{case}_", case, cfg.get("stations_z"), levels=cfg.get("along_s_z"))
    return res


if __name__ == "__main__":
    matplotlib.use("Agg")   # headless when run as a script; imported, the caller keeps its backend
    if len(sys.argv) < 2:
        sys.exit("usage: python hex_post.py <subfolder> [<subfolder> ...]")
    for c in sys.argv[1:]:
        run_case(c)
