#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Diagnostics for sandbox/U_hex_shell.

Streams a one-row-per-step summary to ``output/history.csv``: temperature and
stress at the sample stations of ``case_params.STATIONS``, each averaged over
the six flats, with stresses rotated into the local frame (nn, ss, zz).
``non-regression.py`` reads this CSV.

``__main__`` loads this module automatically when present in the case directory
and calls ``per_step(problem, step, t)`` after every converged step.

Stress is interpolated into DG0, the same cell field the writer exports. Each
station is located on the owned cells of every rank; values are summed with
their hit counts across ranks, so a point on a partition boundary is counted
once per rank that holds it and averaged. The reductions are collective and
run on every rank before the rank-0 write. The file is truncated on the first
call of a run.
"""

import os

import dolfinx
import numpy as np
from mpi4py import MPI

from case_params import NT, STATIONS, flat_axes, flat_point

_CSV = os.path.join(os.path.dirname(__file__), "output", "history.csv")
_NAMES = list(STATIONS)
_NODES = [n for n in _NAMES if STATIONS[n][0] == "node"]
_CELLS = [n for n in _NAMES if STATIONS[n][0] == "cell"]
_HEADER = ",".join(
    ["step", "time_s"]
    + [f"{n}_K" for n in _NODES]
    + [f"{c}_{n}_MPa" for n in _CELLS for c in ("ss", "zz", "nn")]
    + ["ss_flat_spread"]
) + "\n"

_run_started = False
_setup = None


def _locate(problem):
    """Station points (6 per station) and the owned cell holding each one."""
    mesh = problem.mesh
    tdim = mesh.topology.dim
    n_owned = mesh.topology.index_map(tdim).size_local
    tree = dolfinx.geometry.bb_tree(mesh, tdim, entities=np.arange(n_owned, dtype=np.int32))

    pts = np.array([flat_point(m, *STATIONS[n][1:]) for n in _NAMES for m in range(6)])
    hits = dolfinx.geometry.compute_colliding_cells(
        mesh, dolfinx.geometry.compute_collisions_points(tree, pts), pts
    )
    found = np.array([i for i in range(len(pts)) if len(hits.links(i))], dtype=np.int32)
    cells = np.array([hits.links(i)[0] for i in found], dtype=np.int32)

    V_sig = dolfinx.fem.functionspace(mesh, ("DG", 0, (3, 3)))
    return {"pts": pts, "found": found, "cells": cells,
            "sigma": dolfinx.fem.Function(V_sig, name="Stress")}


def _global(fn, comm, n_comp):
    """Value of ``fn`` at every station point, gathered over the ranks."""
    s = _setup
    total = np.zeros((len(s["pts"]), n_comp))
    count = np.zeros(len(s["pts"]))
    if len(s["found"]):
        total[s["found"]] = fn.eval(s["pts"][s["found"]], s["cells"]).reshape(-1, n_comp)
        count[s["found"]] = 1.0
    total = comm.allreduce(total, op=MPI.SUM)
    count = comm.allreduce(count, op=MPI.SUM)
    with np.errstate(invalid="ignore"):
        return total / count[:, None]  # NaN where no rank holds the point


def per_step(problem, step, t):

    global _setup, _run_started
    if _setup is None:
        _setup = _locate(problem)
    comm = problem.mesh.comm

    sigma = next(iter(problem.stress.values()))
    _setup["sigma"].interpolate(dolfinx.fem.Expression(
        sigma, _setup["sigma"].function_space.element.interpolation_points
    ))
    T = _global(problem.T, comm, 1).reshape(len(_NAMES), 6)
    S = _global(_setup["sigma"], comm, 9).reshape(len(_NAMES), 6, 3, 3)

    row = [f"{step}", f"{t:.6e}"]
    for n in _NODES:
        row.append(f"{T[_NAMES.index(n)].mean():.6f}")
    ss_mid = None
    for n in _CELLS:
        Sn = S[_NAMES.index(n)]
        comps = {"ss": [], "zz": [], "nn": []}
        for m in range(6):
            nrm, tang = flat_axes(m)
            comps["ss"].append(tang @ Sn[m] @ tang)
            comps["zz"].append(Sn[m][2, 2])
            comps["nn"].append(nrm @ Sn[m] @ nrm)
        row += [f"{np.mean(comps[c]) / 1e6:.6f}" for c in ("ss", "zz", "nn")]
        if n == f"L{NT - 1}":  # outer layer, mid-flat, mid-height
            ss_mid = np.array(comps["ss"])
    # D6 symmetry: the six flats carry the same stress
    row.append(f"{np.ptp(ss_mid) / abs(ss_mid.mean()):.6e}")

    if comm.rank != 0:
        return

    with open(_CSV, "a" if _run_started else "w") as f:
        if not _run_started:
            f.write(_HEADER)
            _run_started = True
        f.write(",".join(row) + "\n")
