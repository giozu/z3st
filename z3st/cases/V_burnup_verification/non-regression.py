#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: V_burnup_verification

Verifies the burnup state bus (spine.update_state) AND the radial-power source
bus (radial_profile -> set_power f(r,bu)) in a single axisymmetric pellet run.

A solid pellet is held at a constant linear heat rate with a rim-peaking radial
form factor f(r) = 1 + A (r/Ro)^p. The deposited volumetric power is
q(r) = q_avg * f_norm(r), with q_avg = lhr / area, area = pi*Ro^2, and f_norm the
profile normalised to mean 1 (set_power preserves the average rating). Burnup
accumulates as

    bu(r) = q(r) * t / (rho * HM * 8.64e10)   [MWd/kgU],

8.64e10 = 86400 s/day * 1e6 W/MW. Two closed-form checks, independent of each
other:

  1. accumulation magnitude — the nodal-mean burnup equals the flat closed form
         bu_mean = q_avg * t_total / (rho * HM * 8.64e10),
     because the nodal mean of f_norm is 1 by construction. This checks the
     accumulation arithmetic, the unit conversion, and power preservation.
  2. radial shape — the rim/core burnup ratio equals 1 + A = f(Ro)/f(0),
     the *un-normalised* peak factor, independent of the normalisation. This
     checks that the source bus applies the radial shape correctly.

Three figures are produced: accumulation vs time, the radial profile (Z3ST nodal
vs analytical), and a pyvista render of the (r, z) burnup field.
"""

import os
import glob
import yaml
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_verification import *
from z3st.utils.utils_load import generate_power_history
from z3st.materials.fuel_profiles import rim_peaking

CASE_DIR = os.path.dirname(__file__)
OUT = os.path.join(CASE_DIR, "output")
OUT_JSON = os.path.join(OUT, "non-regression.json")

# Burnup is a state that accumulates over the run, so the verification reads the
# *final* step. A multi-step run writes per-step files (fields_NNNN.vtu); a
# single-step run writes fields.vtu — handle both.
_single = os.path.join(OUT, "fields.vtu")
_steps = sorted(glob.glob(os.path.join(OUT, "fields_*.vtu")))
VTU_FILE = _steps[-1] if _steps else _single

# --. geometry + material + history from the case YAML files --..
with open(os.path.join(CASE_DIR, "geometry.yaml")) as f:
    geom = yaml.safe_load(f)
Ro = float(geom["Ro"])
area = np.pi * Ro**2                          # solid pellet cross-section

with open(os.path.join(CASE_DIR, "input.yaml")) as f:
    inp = yaml.safe_load(f)
lhr = float(inp["lhr"][-1])                   # constant linear heat rate (W/m)
t_total = float(max(inp["time"]))             # total irradiation time (s)

with open(os.path.join(CASE_DIR, next(iter(inp["materials"].values())))) as f:
    mat = yaml.safe_load(f)
rho = float(mat["rho"])
hm = float(mat.get("heavy_metal_fraction", 0.8815))
A = float(mat.get("radial_peak_amplitude", 3.0))

# --. analytical references --..
SECONDS_PER_MWD = 8.64e10                      # 86400 s/day * 1e6 W/MW
q_avg = lhr / area                             # mean volumetric power (W/m^3)
BU_MEAN_REF = q_avg * t_total / (rho * hm * SECONDS_PER_MWD)   # area/nodal-mean burnup
RATIO_REF = 1.0 + A                            # rim/core = f(Ro)/f(0)
TOLERANCE = 1e-2

# --. numerical results --..
x, y, z, bu = extract_field(VTU_FILE, field_name="Burnup")
r = np.asarray(x)                              # axisymmetric mesh: x[0] = r
bu = np.asarray(bu)

bu_mean = float(np.mean(bu))
bu_core = float(bu[np.argmin(r)])              # r -> 0
bu_rim = float(bu[np.argmax(r)])               # r = Ro
ratio = bu_rim / max(bu_core, 1e-12)

print(f"[INFO] q_avg = lhr/area = {q_avg:.4e} W/m^3 over t = {t_total:.3e} s")
print(f"[INFO] nodal-mean burnup: numerical = {bu_mean:.6e}, "
      f"analytical = {BU_MEAN_REF:.6e} MWd/kgU")
print(f"[INFO] burnup core = {bu_core:.2f}, rim = {bu_rim:.2f} MWd/kgU")
print(f"[INFO] rim/core ratio: numerical = {ratio:.4f}, analytical (1+A) = {RATIO_REF:.4f}")

errors = {
    "burnup_mean_closed_form": {
        "numerical": bu_mean,
        "reference": BU_MEAN_REF,
        "abs_error": float(abs(bu_mean - BU_MEAN_REF)),
        "rel_error": float(abs(bu_mean - BU_MEAN_REF) / BU_MEAN_REF),
    },
    "rim_core_ratio": {
        "numerical": ratio,
        "reference": RATIO_REF,
        "abs_error": float(abs(ratio - RATIO_REF)),
        "rel_error": float(abs(ratio - RATIO_REF) / RATIO_REF),
    },
}

# --. figure 1: accumulation vs time (per-step mean burnup vs closed form) --..
try:
    times, _, _ = generate_power_history(
        inp["time"], inp["lhr"], n_steps=int(inp["n_steps"]) - 1, filename=None
    )
    bu_step_mean = np.array(
        [float(np.mean(extract_field(f, field_name="Burnup")[3])) for f in _steps]
    )
    if len(bu_step_mean) == len(times):
        t_line = np.linspace(0.0, t_total, 100)
        bu_line = q_avg * t_line / (rho * hm * SECONDS_PER_MWD)
        plt.figure(figsize=(7, 5))
        plt.plot(t_line / 86400.0, bu_line, "k-", lw=2.5, alpha=0.7,
                 label="Analytical  q_avg·t / (ρ·HM·8.64e10)")
        plt.plot(times / 86400.0, bu_step_mean, "o", ms=9, mfc="none", mec="r",
                 mew=2, label="Z3ST mean burnup (update_state)")
        plt.xlabel("time (days)")
        plt.ylabel("mean burnup (MWd/kgU)")
        plt.title("Burnup accumulation: Z3ST vs closed form (constant 20 kW/m)")
        plt.grid(True, ls=":", alpha=0.6)
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(OUT, "burnup_accumulation.png"), dpi=150)
        print("[INFO] burnup_accumulation.png saved")
except Exception as e:
    print(f"[WARNING] accumulation plot skipped: {type(e).__name__}: {e}")

# --. figure 2: radial profile (Z3ST nodal vs analytical) --..
try:
    coords = np.column_stack([r, y, z])
    f_raw = rim_peaking(coords, np.zeros_like(r), mat, model=None)
    f_norm = f_raw / f_raw.mean()             # set_power normalises to mean 1
    bu_profile = q_avg * f_norm * t_total / (rho * hm * SECONDS_PER_MWD)
    order = np.argsort(r)
    plt.figure(figsize=(7, 5))
    plt.plot(r[order] * 1e3, bu_profile[order], "k-", lw=2.5, alpha=0.7,
             label="Analytical  q_avg·f_norm(r)·t / (ρ·HM·8.64e10)")
    plt.scatter(r * 1e3, bu, s=14, facecolors="none", edgecolors="r",
                label="Z3ST nodal burnup")
    plt.axhline(BU_MEAN_REF, color="0.5", ls=":", lw=1.2,
                label=f"mean (flat) = {BU_MEAN_REF:.1f} MWd/kgU")
    plt.xlabel("radius r (mm)")
    plt.ylabel("burnup (MWd/kgU)")
    plt.title("Radial burnup distribution: rim-peaking form factor")
    plt.grid(True, ls=":", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "burnup_radial_profile.png"), dpi=150)
    print("[INFO] burnup_radial_profile.png saved")
except Exception as e:
    print(f"[WARNING] radial-profile plot skipped: {type(e).__name__}: {e}")

# --. figure 3: pyvista field render (off-screen) --..
try:
    import pyvista as pv
    pv.OFF_SCREEN = True
    grid = pv.read(VTU_FILE)
    p = pv.Plotter(off_screen=True, window_size=(720, 900))
    sbar = {
        "title": "Burnup (MWd/kgU)", "vertical": True,
        "position_x": 0.80, "position_y": 0.15, "width": 0.12, "height": 0.7,
        "title_font_size": 18, "label_font_size": 14,
    }
    p.add_mesh(grid, scalars="Burnup", cmap="inferno", show_edges=False,
               scalar_bar_args=sbar)
    p.view_xy()
    p.enable_parallel_projection()
    p.add_text("axis  r ->  surface", position="lower_edge", font_size=9)
    p.camera.zoom(1.2)
    p.screenshot(os.path.join(OUT, "burnup_field_pyvista.png"))
    print("[INFO] burnup_field_pyvista.png saved")
except Exception as e:
    print(f"[WARNING] pyvista render skipped: {type(e).__name__}: {e}")

pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] non-regression completed.\n")
