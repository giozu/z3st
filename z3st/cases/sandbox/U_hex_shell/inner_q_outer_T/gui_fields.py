#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Grid loader for ``interactive_gui.ipynb``, reading only ``output/``.

D, H and the hexagon orientation come from the mesh stored in
``output/fields.h5``, not from the case YAML files, so the notebook shows a
past run whatever the inputs currently say.
"""

import os

import h5py
import numpy as np
import pyvista as pv

from z3st.utils.utils_extract_xdmf import extract_field_xdmf

HERE = os.path.dirname(os.path.abspath(__file__))
XDMF = os.path.join(HERE, "output", "fields.xdmf")

with h5py.File(XDMF.replace(".xdmf", ".h5"), "r") as _f:
    _PTS = np.array(_f["Mesh/mesh/geometry"])
    _TOPO = np.array(_f["Mesh/mesh/topology"])

_r = np.hypot(_PTS[:, 0], _PTS[:, 1])
_corner = _PTS[np.argmax(_r)]
D = 2 * _r.max()                                    # (m) outer circumscribed diameter
H = _PTS[:, 2].max() - _PTS[:, 2].min()             # (m) height
A0 = np.arctan2(_corner[1], _corner[0]) % (np.pi / 3)  # corner angle, mod 60 deg


def flat_frame(x, y):
    """Outward normal and tangent of the flat each point lies on."""
    sector = np.floor((np.arctan2(y, x) - A0) / (np.pi / 3))
    phi_n = A0 + sector * np.pi / 3 + np.pi / 6
    zero = np.zeros_like(phi_n)
    n = np.stack([np.cos(phi_n), np.sin(phi_n), zero], axis=1)
    tang = np.stack([-np.sin(phi_n), np.cos(phi_n), zero], axis=1)
    return n, tang


def load_grid(step_index=-1):
    """pyvista grid of a step: T, u, |u| on the nodes; stress and strain
    components in the local frame of the flats (nn, ss, zz) and the von Mises
    stress on the cells."""
    grid = pv.UnstructuredGrid(
        np.hstack([np.full((len(_TOPO), 1), _TOPO.shape[1]), _TOPO]).ravel(),
        np.full(len(_TOPO), pv.CellType.HEXAHEDRON, dtype=np.uint8), _PTS)

    get = lambda name: extract_field_xdmf(XDMF, name, step_index, return_coords=False)
    U = get("Displacement")
    grid.point_data["T (K)"] = get("Temperature")
    grid.point_data["u"] = U
    grid.point_data["u (um)"] = np.linalg.norm(U, axis=1) * 1e6

    centres = _PTS[_TOPO].mean(axis=1)
    n, tang = flat_frame(centres[:, 0], centres[:, 1])
    ez = np.broadcast_to([0.0, 0.0, 1.0], n.shape)
    for name, sym, scale, unit in (("Stress", "sigma", 1e-6, "MPa"),
                                   ("Strain", "eps", 1.0, "-")):
        A = get(name).reshape(-1, 3, 3)
        for comp, e in (("nn", n), ("ss", tang), ("zz", ez)):
            grid.cell_data[f"{sym}_{comp} ({unit})"] = (
                np.einsum("ci,cij,cj->c", e, A, e) * scale)
        if sym == "sigma":
            dev = A - np.trace(A, axis1=1, axis2=2)[:, None, None] / 3 * np.eye(3)
            grid.cell_data["sigma_vm (MPa)"] = (
                np.sqrt(1.5 * np.einsum("cij,cij->c", dev, dev)) * scale)
    return grid


def cells_as_points(mesh):
    """Outer surface with every cell field copied onto unshared nodes.

    Mesa llvmpipe (WSL, also under ParaView) draws nothing for cell scalars
    above a few hundred cells; point scalars on per-cell nodes render the same
    piecewise-constant field.
    """
    s = mesh.extract_surface().cast_to_unstructured_grid().separate_cells()
    owner = np.empty(s.n_points, dtype=int)
    owner[s.cell_connectivity] = np.repeat(np.arange(s.n_cells), np.diff(s.cell_offsets))
    for name in s.cell_data.keys():
        s.point_data[name] = s.cell_data[name][owner]
    return s
