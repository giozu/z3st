#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Diagnostics for regression/fg_test_fuel.

Streams a per-step summary to ``output/history.csv`` (burnup, temperatures,
fission-gas scalars) for ``plots.py``. ``__main__`` calls ``per_step`` after
every converged step.
"""

import os
import numpy as np

_CSV = os.path.join(os.path.dirname(__file__), "output", "history.csv")
_HEADER = ("step,time_s,time_days,burnup_avg_MWdkgU,burnup_max_MWdkgU,"
           "T_max_K,T_min_K,"
           "fg_produced_avg,fg_ingrain_avg,fg_gb_avg,fg_released_avg,fgr_frac\n")


def _fg_averages(problem):
    """Fuel-average total fission gas (Xe+Kr) [at/m^3] and fractional release,
    over the fissile dofs (``problem._sciantix_dofs``). Zeros when coupling off."""
    fields = getattr(problem, "fg_fields", None)
    dofs = getattr(problem, "_sciantix_dofs", None)
    if not fields or dofs is None or len(dofs) == 0:
        return 0.0, 0.0, 0.0, 0.0, 0.0
    prod = np.asarray(fields["produced"].x.array)[dofs]
    avgs = [float(np.asarray(fields[k].x.array)[dofs].mean())
            for k in ("produced", "in_grain", "grain_boundary", "released")]
    rel_sum = float(np.asarray(fields["released"].x.array)[dofs].sum())
    fgr = rel_sum / float(prod.sum()) if prod.sum() > 0 else 0.0
    return (*avgs, fgr)


def per_step(problem, step, t):

    bu_avg = bu_max = 0.0
    bu_fn = getattr(problem, "burnup", None)
    if bu_fn is not None:
        bu = np.asarray(bu_fn.x.array)
        if np.any(bu > 0):
            bu_avg = float(bu[bu > 0].mean())
            bu_max = float(bu.max())

    T_fn = getattr(problem, "T", None)
    t_max = float(np.asarray(T_fn.x.array).max()) if T_fn is not None else float("nan")
    t_min = float(np.asarray(T_fn.x.array).min()) if T_fn is not None else float("nan")

    fg_prod, fg_ingrain, fg_gb, fg_rel, fgr = _fg_averages(problem)

    # truncate at step 0 so a stale CSV can't graft a second burnup sweep onto this one
    fresh = (step == 0)
    with open(_CSV, "w" if fresh else "a") as f:
        if fresh or not os.path.exists(_CSV):
            f.write(_HEADER)
        f.write(f"{step},{t:.6e},{t / 86400.0:.4f},{bu_avg:.6e},{bu_max:.6e},"
                f"{t_max:.2f},{t_min:.2f},"
                f"{fg_prod:.6e},{fg_ingrain:.6e},{fg_gb:.6e},{fg_rel:.6e},{fgr:.6f}\n")
