#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.3.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Diagnostics for verification/fuel/creep_shrink_fit_2D.

Streams a one-row-per-step summary to ``output/history.csv`` of scalars the
solver has already computed. ``plots.py`` reads this CSV for the PCMI curves.

``__main__`` loads this module automatically when present in the case directory
and calls ``per_step(problem, step, t)`` after every converged step.

Under MPI the field statistics are reduced across ranks and only rank 0 writes,
so the CSV holds one row per step with global values. The reductions are
collective and must therefore run on every rank before the rank-0 write. The
file is truncated on the first call of a run rather than appended to.
"""

import os
import numpy as np
from mpi4py import MPI

from case_params import check_consistency

# At import, i.e. before step 0. SystemExit is not an Exception, so __main__'s
# loader does not downgrade it to a warning.
if _problems := check_consistency():
    raise SystemExit("[diagnostics] configuration is incoherent with Esposito "
                     "eq. (21):\n  - " + "\n  - ".join(_problems))

_CSV = os.path.join(os.path.dirname(__file__), "output", "history.csv")
_HEADER = ("step,time_s,time_days,burnup_avg_MWdkgU,burnup_max_MWdkgU,"
           "gap_um,contact_pressure_MPa,T_max_K,T_min_K\n")

_run_started = False


def _owned(fn):
    """Local dof values excluding ghosts, so a reduction over ranks counts
    every entry exactly once."""
    imap = fn.function_space.dofmap.index_map
    bs = fn.function_space.dofmap.index_map_bs
    return np.asarray(fn.x.array)[: imap.size_local * bs]


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

    if comm is not None and comm.rank != 0:
        return

    global _run_started
    with open(_CSV, "a" if _run_started else "w") as f:
        if not _run_started:
            f.write(_HEADER)
            _run_started = True
        f.write(f"{step},{t:.6e},{t / 86400.0:.4f},{bu_avg:.6e},{bu_max:.6e},"
                f"{gap_um:.4f},{p_mpa:.4f},{t_max:.2f},{t_min:.2f}\n")
