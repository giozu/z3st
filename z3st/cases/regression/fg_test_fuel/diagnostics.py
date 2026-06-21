#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Diagnostics for regression/fg_test_fuel.

Streams a per-step summary to ``output/history.csv`` (burnup, temperatures,
fission-gas scalars) for ``plots.py``. ``__main__`` calls ``per_step`` on every
rank after each converged step; the reductions are collective (global) and only
rank 0 writes.
"""

import os
import numpy as np
from mpi4py import MPI

_CSV = os.path.join(os.path.dirname(__file__), "output", "history.csv")
_HEADER = ("step,time_s,time_days,burnup_avg_MWdkgU,burnup_max_MWdkgU,"
           "T_max_K,T_min_K,"
           "fg_produced_avg,fg_ingrain_avg,fg_gb_avg,fg_released_avg,fgr_frac\n")


def _owned(fn):
    """Number of owned (non-ghost) dofs of a scalar Function."""
    return fn.function_space.dofmap.index_map.size_local


def _fg_averages(problem):
    """Global fuel-average total fission gas (Xe+Kr) (at/m^3) and fractional
    release over the fissile dofs. Collective; zeros when the coupling is off."""
    comm = problem.mesh.comm
    fields = getattr(problem, "fg_fields", None)
    dofs = getattr(problem, "_sciantix_dofs", None)
    if not fields or dofs is None:   # coupling on/off is global -> no rank mismatch
        return 0.0, 0.0, 0.0, 0.0, 0.0
    dofs = np.asarray(dofs)
    dofs_own = dofs[dofs < _owned(fields["produced"])]   # owned only: ghosts would double-count
    cnt = comm.allreduce(int(dofs_own.size), op=MPI.SUM)

    def gsum(key):
        return comm.allreduce(float(np.asarray(fields[key].x.array)[dofs_own].sum()), op=MPI.SUM)

    sums = {k: gsum(k) for k in ("produced", "in_grain", "grain_boundary", "released")}
    avgs = [sums[k] / cnt if cnt > 0 else 0.0
            for k in ("produced", "in_grain", "grain_boundary", "released")]
    fgr = sums["released"] / sums["produced"] if sums["produced"] > 0 else 0.0
    return (*avgs, fgr)


def per_step(problem, step, t):
    comm = problem.mesh.comm

    bu_avg = bu_max = 0.0
    bu_fn = getattr(problem, "burnup", None)
    if bu_fn is not None:
        bu = np.asarray(bu_fn.x.array)[:_owned(bu_fn)]
        pos = bu > 0
        s = comm.allreduce(float(bu[pos].sum()), op=MPI.SUM)
        c = comm.allreduce(int(pos.sum()), op=MPI.SUM)
        bu_avg = s / c if c > 0 else 0.0
        bu_max = comm.allreduce(float(bu.max()) if bu.size else 0.0, op=MPI.MAX)

    T_fn = getattr(problem, "T", None)
    if T_fn is not None:
        Ta = np.asarray(T_fn.x.array)   # owned+ghost ok for max/min
        t_max = comm.allreduce(float(Ta.max()), op=MPI.MAX)
        t_min = comm.allreduce(float(Ta.min()), op=MPI.MIN)
    else:
        t_max = t_min = float("nan")

    fg_prod, fg_ingrain, fg_gb, fg_rel, fgr = _fg_averages(problem)

    if comm.rank != 0:   # all ranks reduced the same globals; only rank 0 writes
        return
    fresh = (step == 0)   # truncate at step 0 so a stale CSV can't graft a second sweep
    with open(_CSV, "w" if fresh else "a") as f:
        if fresh or not os.path.exists(_CSV):
            f.write(_HEADER)
        f.write(f"{step},{t:.6e},{t / 86400.0:.4f},{bu_avg:.6e},{bu_max:.6e},"
                f"{t_max:.2f},{t_min:.2f},"
                f"{fg_prod:.6e},{fg_ingrain:.6e},{fg_gb:.6e},{fg_rel:.6e},{fgr:.6f}\n")
