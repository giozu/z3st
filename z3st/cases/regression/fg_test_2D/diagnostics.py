#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Diagnostics for regression/fg_test_2D.

Streams a one-row-per-step summary to ``output/history.csv`` of scalars the
solver has already computed (burnup, gap and contact pressure, temperatures,
fission-gas concentrations). ``plots.py`` reads this CSV; field snapshots come
from VTU (short runs) or ParaView (the single XDMF).

``__main__`` loads this module automatically when present in the case directory
and calls ``per_step(problem, step, t)`` after every converged step.

Under MPI the field statistics are reduced across ranks and only rank 0 writes,
so the CSV holds one row per step with global values. The reductions are
collective and must therefore run on every rank before the rank-0 write. The
file is truncated on the first call of a run rather than appended to, so a
re-run without ``Allclean`` replaces the previous history instead of
concatenating with it.
"""

import os
import numpy as np
from mpi4py import MPI

_CSV = os.path.join(os.path.dirname(__file__), "output", "history.csv")
_HEADER = ("step,time_s,time_days,burnup_avg_MWdkgU,burnup_max_MWdkgU,"
           "gap_um,contact_pressure_MPa,T_max_K,T_min_K,"
           "fg_produced_avg,fg_ingrain_avg,fg_gb_avg,fg_released_avg,fgr_frac\n")

_run_started = False


def _n_owned(fn):
    """Number of owned (non-ghost) dof entries of a Function."""
    dofmap = fn.function_space.dofmap
    return dofmap.index_map.size_local * dofmap.index_map_bs


def _owned(fn):
    """Local dof values excluding ghosts, so a reduction over ranks counts
    every entry exactly once."""
    return np.asarray(fn.x.array)[: _n_owned(fn)]


def _fg_averages(problem, comm, parallel):
    """Global fuel-average total fission-gas (Xe+Kr) concentrations [at/m^3] and
    the fuel-average fractional release. Returns zeros when the coupling is off.

    The FG Functions are nonzero only on the fissile dofs, so the average is taken
    over those dofs (``problem._sciantix_dofs``); the fractional release uses the
    summed produced/released over the fuel (a volume-consistent average since all
    fissile dofs share the V_t element).
    """
    fields = getattr(problem, "fg_fields", None)
    dofs = getattr(problem, "_sciantix_dofs", None)
    if not fields or dofs is None:   # coupling on/off is global -> no rank mismatch
        return 0.0, 0.0, 0.0, 0.0, 0.0

    dofs = np.asarray(dofs)
    dofs_own = dofs[dofs < _n_owned(fields["produced"])]   # ghosts would double-count
    keys = ("produced", "in_grain", "grain_boundary", "released")

    cnt = int(dofs_own.size)
    sums = {k: float(np.asarray(fields[k].x.array)[dofs_own].sum()) for k in keys}
    if parallel:
        cnt = comm.allreduce(cnt, op=MPI.SUM)
        for k in keys:
            sums[k] = comm.allreduce(sums[k], op=MPI.SUM)

    avgs = [sums[k] / cnt if cnt > 0 else 0.0 for k in keys]
    fgr = sums["released"] / sums["produced"] if sums["produced"] > 0 else 0.0
    return (*avgs, fgr)


def per_step(problem, step, t):

    comm = getattr(getattr(problem, "mesh", None), "comm", None)
    parallel = comm is not None and comm.size > 1

    bu_avg = bu_max = 0.0
    bu_fn = getattr(problem, "burnup", None)
    if bu_fn is not None:
        bu = _owned(bu_fn)
        pos = bu[bu > 0]
        bu_sum, bu_cnt = float(pos.sum()), int(pos.size)
        bu_peak = float(bu.max()) if bu.size else 0.0
        if parallel:
            bu_sum = comm.allreduce(bu_sum, op=MPI.SUM)
            bu_cnt = comm.allreduce(bu_cnt, op=MPI.SUM)
            bu_peak = comm.allreduce(bu_peak, op=MPI.MAX)
        if bu_cnt:
            bu_avg = bu_sum / bu_cnt
            bu_max = bu_peak

    # Both are single global scalars maintained by the contact model, so they
    # need no reduction.
    gap_um = float(getattr(problem, "_last_gap", float("nan")) * 1e6)
    p_mpa = float(getattr(problem, "_last_pressure", 0.0) / 1e6)

    T_fn = getattr(problem, "T", None)
    if T_fn is None:
        t_max = t_min = float("nan")
    else:
        T = _owned(T_fn)
        t_max = float(T.max()) if T.size else -np.inf
        t_min = float(T.min()) if T.size else np.inf
        if parallel:
            t_max = comm.allreduce(t_max, op=MPI.MAX)
            t_min = comm.allreduce(t_min, op=MPI.MIN)

    fg_prod, fg_ingrain, fg_gb, fg_rel, fgr = _fg_averages(problem, comm, parallel)

    if comm is not None and comm.rank != 0:
        return

    global _run_started
    with open(_CSV, "a" if _run_started else "w") as f:
        if not _run_started:
            f.write(_HEADER)
            _run_started = True
        f.write(f"{step},{t:.6e},{t / 86400.0:.4f},{bu_avg:.6e},{bu_max:.6e},"
                f"{gap_um:.4f},{p_mpa:.4f},{t_max:.2f},{t_min:.2f},"
                f"{fg_prod:.6e},{fg_ingrain:.6e},{fg_gb:.6e},{fg_rel:.6e},{fgr:.6f}\n")
