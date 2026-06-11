#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST case-local diagnostics --.. ..- .-.. .-.. ---
"""
PCMI diagnostics for regression/pwr_rod_2D.

Streams a one-row-per-step summary to ``output/history.csv``. It records scalars
the solver has already computed — the contact model's last gap/pressure
(``problem._last_gap`` (m), ``problem._last_pressure`` (Pa)), the burnup field,
and the temperature field — so it works whether the run writes VTU or XDMF and
needs no per-step VTU files at all. ``plots.py`` reads this CSV for the PCMI
curves; field snapshots come from VTU (short runs) or ParaView (the single XDMF).

``__main__`` loads this module automatically when present in the case directory
and calls ``per_step(problem, step, t)`` after every converged step.
"""

import os
import numpy as np

_CSV = os.path.join(os.path.dirname(__file__), "output", "history.csv")
_HEADER = ("step,time_s,time_days,burnup_avg_MWdkgU,burnup_max_MWdkgU,"
           "gap_um,contact_pressure_MPa,T_max_K\n")


def per_step(problem, step, t):
    # Fuel-average and peak burnup (the field is nonzero only in the fissile
    # pellet, so the >0 mask selects the fuel region).
    bu_avg = bu_max = 0.0
    bu_fn = getattr(problem, "burnup", None)
    if bu_fn is not None:
        bu = np.asarray(bu_fn.x.array)
        if np.any(bu > 0):
            bu_avg = float(bu[bu > 0].mean())
            bu_max = float(bu.max())

    # Contact gap / pressure as stored by the contact model.
    gap_um = float(getattr(problem, "_last_gap", float("nan")) * 1e6)
    p_mpa = float(getattr(problem, "_last_pressure", 0.0) / 1e6)

    # Peak temperature.
    T_fn = getattr(problem, "T", None)
    t_max = float(np.asarray(T_fn.x.array).max()) if T_fn is not None else float("nan")

    is_new = not os.path.exists(_CSV)
    with open(_CSV, "a") as f:
        if is_new:
            f.write(_HEADER)
        f.write(f"{step},{t:.6e},{t / 86400.0:.4f},{bu_avg:.6e},{bu_max:.6e},"
                f"{gap_um:.4f},{p_mpa:.4f},{t_max:.2f}\n")
