#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
# Author: Bianca Funaro

"""
Z3ST case: sandbox/U_hex_shell

Hexagonal steel shell (assembly wrapper) with a through-wall temperature
gradient: Neumann heat flux q on the inner faces, Dirichlet T_o on the outer
faces, no mechanical load (rigid-body modes removed by the null space); an
axial clamp at the bottom (u_z = 0) is a mirror plane and stays within the
validity range. The derivation is in
"Exercise - Thermal stress in flat plate.ipynb".

All metrics are read from ``output/history.csv`` — the per-step values at the
sample stations of ``case_params.STATIONS``, streamed by the case-local
``diagnostics.py`` and averaged over the six flats — so the check is
independent of the output format.

Away from the corners and the ends, each flat is a slab with a linear
temperature profile. A free plate would take it entirely into bending, with
zero stress; the closed section suppresses the bending, and the stresses are

    k((T_i + T_o)/2) (T_i - T_o) = q t                         (linear k, exact)
    sigma_nn = 0
    sigma_ss = sigma_zz = alpha E / (1 - nu) (T_mean - T(xi))  (odd in xi)

Metrics with a closed form, along the three lines of the local frame:

  * xi (through the wall), mid-flat, mid-height: the wall drop, sigma_ss and
    sigma_zz at the nt cell layers, sigma_nn = 0, and the symmetries of the
    solution: sigma_ss odd in xi, equal on the six flats (D6).
    sigma_zz carries a uniform membrane offset: axial force balance is over
    the whole section, corners included, which the slab model does not see.
    The offset is tracked, and sigma_zz is checked with it removed.
  * s (along the flat), mid-height: sigma_ss uniform up to |s| = L/4.
  * z (along the axis), mid-flat: T_i uniform over the height, sigma_ss
    uniform up to a quarter of the height from the end.

The corner and end disturbances have no closed form — they are recorded with
``tracked`` so the analytic pass/fail gate ignores them, and are protected
purely by the gold comparison (``regression_check`` vs
``output/non-regression_gold.json``).

The thermal problem is solved as a repeated steady state with k(T) taken from
the previous step, so the wall drop of the last two steps must agree.
"""

import os
import csv
import yaml
import numpy as np

from case_params import DT, NT, SIGMA_SURF, T_I, XI_CELLS, sigma_wall
from z3st.utils.non_regression import error_metric, finish, metric, tracked
from z3st.utils.utils_load import generate_power_history

CASE_DIR = os.path.dirname(__file__)
OUT = os.path.join(CASE_DIR, "output")
OUT_JSON = os.path.join(OUT, "non-regression.json")
HISTORY = os.path.join(OUT, "history.csv")

TOLERANCE = 1e-2

# --. trajectory from the case-local diagnostics CSV --..
with open(HISTORY) as f:
    rows = list(csv.DictReader(f))
if not rows:
    raise RuntimeError(f"{HISTORY} is empty — did the run complete?")
last = rows[-1]
if any(np.isnan(float(v)) for v in last.values()):
    raise RuntimeError("history.csv holds NaN: a station lies outside the mesh — "
                       "geometry.yaml does not match the solved mesh")

# --. expected number of steps from the input time grid --..
with open(os.path.join(CASE_DIR, "input.yaml")) as f:
    inp = yaml.safe_load(f)
raw_n_steps = inp["n_steps"]
n_increments = raw_n_steps if isinstance(raw_n_steps, (list, tuple)) else int(raw_n_steps) - 1
times, _, _ = generate_power_history(inp["time"], inp["lhr"], n_steps=n_increments, filename=None)
if len(rows) != len(times):
    print(f"[WARNING] history.csv has {len(rows)} rows, the time grid {len(times)} steps.")

# --. line 1: through the wall (xi) --..
dT = float(last["T_i_K"]) - float(last["T_o_K"])
dT_prev = float(rows[-2]["T_i_K"]) - float(rows[-2]["T_o_K"]) if len(rows) > 1 else dT
ss = np.array([float(last[f"ss_L{j}_MPa"]) for j in range(NT)]) * 1e6
zz = np.array([float(last[f"zz_L{j}_MPa"]) for j in range(NT)]) * 1e6
nn = np.array([float(last[f"nn_L{j}_MPa"]) for j in range(NT)]) * 1e6
ref = sigma_wall(XI_CELLS)
zz_offset = zz.mean()
flat_spread = float(last["ss_flat_spread"])

print(f"[INFO] wall drop T_i - T_o : numerical = {dT:.5f}, analytic = {DT:.5f} K")
for j in range(NT):
    print(f"[INFO] xi = {XI_CELLS[j] * 1e3:+.3f} mm : analytic = {ref[j] / 1e6:+.4f}, "
          f"ss = {ss[j] / 1e6:+.4f}, zz = {zz[j] / 1e6:+.4f}, nn = {nn[j] / 1e6:+.4f} MPa")
print(f"[INFO] sigma_zz membrane offset : {zz_offset / 1e6:+.4f} MPa")
print(f"[INFO] flat-to-flat spread of sigma_ss (D6) : {flat_spread:.2e}")

# --. line 2: along the flat (s) --..
sigma_outer = sigma_wall(XI_CELLS[-1])
ss_window_s = float(last["ss_window_s_MPa"]) * 1e6
corner_ratio = float(last["ss_corner_MPa"]) * 1e6 / sigma_outer
corner_drop = float(last["T_i_K"]) - float(last["T_i_corner_K"])

print(f"[INFO] s line : sigma_ss(L/4) = {ss_window_s / 1e6:+.4f} MPa (analytic "
      f"{sigma_outer / 1e6:+.4f}), corner / analytic = {corner_ratio:.3f}, "
      f"corner T_i drop = {corner_drop:.4f} K")

# --. line 3: along the axis (z) --..
ss_window_z = float(last["ss_window_z_MPa"]) * 1e6
bottom_ratio = float(last["ss_bottom_MPa"]) * 1e6 / sigma_outer
top_ratio = float(last["ss_top_MPa"]) * 1e6 / sigma_outer
T_i_z = np.array([float(last[k]) for k in ("T_i_bottom_K", "T_i_K", "T_i_top_K")])

print(f"[INFO] z line : sigma_ss(H/4) = {ss_window_z / 1e6:+.4f} MPa, "
      f"bottom / analytic = {bottom_ratio:.3f}, top / analytic = {top_ratio:.3f}")

errors = {
    # xi
    "wall_dT_K": metric(dT, DT),
    "wall_dT_last_step_change": error_metric(abs(dT - dT_prev) / DT),
    "sigma_ss_xi_L2": error_metric(np.linalg.norm(ss - ref) / np.linalg.norm(ref)),
    "sigma_zz_xi_L2": error_metric(np.linalg.norm(zz - zz_offset - ref) / np.linalg.norm(ref)),
    "sigma_nn_xi_max": error_metric(np.abs(nn).max() / SIGMA_SURF),
    "sigma_ss_xi_odd": error_metric(abs(ss[0] + ss[-1]) / abs(ss[-1])),
    "sigma_ss_flat_spread": error_metric(flat_spread),
    "sigma_zz_offset_MPa": tracked(zz_offset / 1e6),
    # s
    "sigma_ss_s_window": metric(ss_window_s, sigma_outer),
    "sigma_ss_corner_ratio": tracked(corner_ratio),
    "T_inner_corner_drop_K": tracked(corner_drop),
    # z
    "T_inner_z_max": error_metric(np.abs(T_i_z - T_i_z[1]).max() / DT),
    "sigma_ss_z_window": metric(ss_window_z, sigma_outer),
    "sigma_ss_bottom_ratio": tracked(bottom_ratio),
    "sigma_ss_top_ratio": tracked(top_ratio),
}

finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)
