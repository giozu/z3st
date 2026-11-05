# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# beam1d: Euler-Bernoulli model of the hexagonal wrapper
# Version: 0.2 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Temperature T(x, y, z) and thermal eigenstrain eps*(x) = alpha(T) (T - T_ref)
from the thermal block of boundary_conditions.yaml (one entry, the whole wall):

    type: expression   temperature: "<expr in x, y, z and np>"  (K)
                       the z3st syntax (models/thermal_model.py), evaluated with
                       eval on numpy arrays: the case input is trusted.
    type: table        file: CSV with a z column and one column per key point
                       z_column, z_unit (m | mm), unit (K | C)
                       points: {corners: [6 names], flats: [6 names]}
                         corner k of beam1d.section and the mid-point of flat k
                       section: contour_linear | ma12

A table gives T at 12 key points of the mid-wall contour (Ma et al. 2019,
Fig. 8: A1-A6 at the corners, S1-S6 at mid-flat), linear in z between the
rows and constant beyond each column's data. Around the section:

    contour_linear   T linear along each half flat, corner -> mid-flat -> corner
    ma12             12 sectors of constant T (Ma Fig. 2): the central half of
                     each flat takes S, the quarter flats on both sides of a
                     corner take A (the widths are not given in the paper: TBC)

Both are uniform through the wall. The section splits each flat at its
quarter and mid points (beam1d.section.SPLITS), so the section integrals of
either model are exact.

alpha in the material card is a number (1/K) or a table {T: [K], value: [1/K]}
of the secant coefficient from T_ref (Ma Eq. 30 uses alpha_ij (T_ij - T_f)),
interpolated linearly in T and constant beyond the table.
"""

import os

import numpy as np


def alpha_function(alpha):
    """alpha(T) (1/K) from the card entry: a number or a secant table."""
    if isinstance(alpha, dict):
        T_tab = np.asarray(alpha["T"], dtype=float)
        a_tab = np.asarray(alpha["value"], dtype=float)
        return lambda T: np.interp(T, T_tab, a_tab)
    a = float(alpha)
    return lambda T: a + 0.0 * T


def expression(expr):
    """T(x) for x = (x, y, z) from a z3st expression string."""
    def T(x):
        X, Y, Z = np.broadcast_arrays(*x)
        return np.broadcast_to(eval(expr, {"__builtins__": {}, "np": np}, {"x": X, "y": Y, "z": Z}), X.shape)
    return T


class KeyPointTable:
    """T at the 12 key points of the contour vs z (see the module docstring).

    columns: key point -> (z (m), T (K)) arrays; from_csv reads them from a file.
    """

    def __init__(self, columns, section, points, model="contour_linear"):
        if model not in ("contour_linear", "ma12"):
            raise ValueError(f"unknown section model '{model}' (contour_linear | ma12)")
        self.section, self.model = section, model
        self.names = {"corners": list(points["corners"]), "flats": list(points["flats"])}
        self.columns = {n: (np.asarray(columns[n][0], dtype=float), np.asarray(columns[n][1], dtype=float))
                        for n in self.names["corners"] + self.names["flats"]}
        self.z_range = (min(c[0].min() for c in self.columns.values()),
                        max(c[0].max() for c in self.columns.values()))

    @classmethod
    def from_csv(cls, path, section, points, z_column="z_mm", z_unit="mm", unit="C", model="contour_linear"):
        data = np.genfromtxt(path, delimiter=",", names=True)
        z_scale = {"m": 1.0, "mm": 1e-3}[z_unit]
        T_shift = {"K": 0.0, "C": 273.15}[unit]
        columns = {}
        for name in list(points["corners"]) + list(points["flats"]):
            ok = np.isfinite(data[name]) & np.isfinite(data[z_column])
            columns[name] = (data[z_column][ok] * z_scale, data[name][ok] + T_shift)
        table = cls(columns, section, points, model)
        table.path = path
        return table

    def at(self, name, z):
        """T (K) of one key point at the heights z."""
        zt, Tt = self.columns[name]
        return np.interp(z, zt, Tt)

    def __call__(self, x):
        X, Y, Z = np.broadcast_arrays(*(np.asarray(a, dtype=float) for a in x))
        k, f = self.section.contour_coords(X, Y)
        Tc = np.stack([self.at(n, Z) for n in self.names["corners"]])   # (6, ...)
        Tf = np.stack([self.at(n, Z) for n in self.names["flats"]])
        take = lambda T, j: np.take_along_axis(T, j[None, ...], axis=0)[0]   # noqa: E731
        T0, Tm, T1 = take(Tc, k), take(Tf, k), take(Tc, (k + 1) % 6)
        if self.model == "contour_linear":
            return np.where(f <= 0.5, T0 + (Tm - T0) * 2.0 * f, Tm + (T1 - Tm) * (2.0 * f - 1.0))
        return np.where(f < 0.25, T0, np.where(f <= 0.75, Tm, T1))


def build(case, section, section_model=None):
    """T(x) (K) and eps*(x) of the case, both callables of x = (x, y, z).

    section_model overrides the 'section' of a table (contour_linear | ma12).
    """
    entries = case.bcs.get("thermal", {}).get(case.mat_name, [])
    if len(entries) != 1:
        raise ValueError(f"thermal: {case.mat_name}: one entry expected (the whole wall), found {len(entries)}")
    bc = entries[0]
    if bc["type"] == "expression":
        T = expression(bc["temperature"])
    elif bc["type"] == "table":
        T = KeyPointTable.from_csv(os.path.join(case.dir, bc["file"]), section, bc["points"],
                                   bc.get("z_column", "z_mm"), bc.get("z_unit", "mm"), bc.get("unit", "C"),
                                   section_model or bc.get("section", "contour_linear"))
    else:
        raise ValueError(f"unknown thermal type '{bc['type']}' (expression | table)")
    alpha = alpha_function(case.mat["alpha"])
    T_ref = float(case.mat["T_ref"])

    def eps(x):
        Tx = T(x)
        return alpha(Tx) * (Tx - T_ref)
    return T, eps
