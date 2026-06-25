#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: verification/fuel/porosity_migration

Reproduces the as-fabricated porosity migration test case of Barani et al.
(2022), J. Nucl. Mater. 558, 153341 (originally Novascone et al. 2018, ref. 14):
a (U,Pu)O2 pellet sector brought to 500 W/cm with the outer rim ramped 623->1300 K
over 1e4 s. Pores migrate up the thermal gradient and pile up at the centre,
forming a central void surrounded by a restructured (columnar) low-porosity zone.

The absolute checks are anchored to quantities the paper actually reports, NOT to
values produced by this run:
  1. A central void forms: p(r=0) ~ 1.0.
  2. The central void extends to ~0.2 of the relative radius (Barani Fig. 4: the
     thermal gradient and pore velocity vanish inside ~0.2 r/Ro because p = 1).
  3. The outer rim keeps its fabrication porosity p ~ 0.15. This is now a genuine
     prediction (no Dirichlet inflow is imposed; the rim stays at 0.15 because the
     pore velocity there is negligible).

The centre temperature is reported for information and tracked by the gold
regression, but it is not asserted against a fabricated reference (the paper does
not tabulate it).
"""

import os
import glob
import yaml
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from z3st.utils.utils_extract_vtu import extract_field
from z3st.utils.utils_extract_xdmf import extract_field_xdmf
from z3st.utils.utils_verification import pass_fail_check, regression_check

CASE_DIR = os.path.dirname(__file__)
OUT = os.path.join(CASE_DIR, "output")
OUT_JSON = os.path.join(OUT, "non-regression.json")

# Output is VTU in serial (Allrun) and XDMF in parallel (Allrun_mpi): the writer
# switches vtu -> xdmf automatically when comm.size > 1. Read whichever exists.
_vtu_steps = sorted(glob.glob(os.path.join(OUT, "fields_*.vtu")))
_vtu_single = os.path.join(OUT, "fields.vtu")
_xdmf = os.path.join(OUT, "fields.xdmf")

if _vtu_steps or os.path.exists(_vtu_single):
    _vtu_file = _vtu_steps[-1] if _vtu_steps else _vtu_single
    def get_field(name):
        return extract_field(_vtu_file, field_name=name)
elif os.path.exists(_xdmf):
    def get_field(name):
        return extract_field_xdmf(_xdmf, field_name=name, step_index=-1)
else:
    raise FileNotFoundError(
        f"No fields.vtu / fields_*.vtu / fields.xdmf found in {OUT}"
    )

# --. Load geometry --..
with open(os.path.join(CASE_DIR, "geometry.yaml")) as f:
    geom = yaml.safe_load(f)
Ro = float(geom["Ro"])

# --. Extract fields --..
x, y, z, p_vals = get_field("Porosity")
_, _, _, T_vals = get_field("Temperature")

x = np.asarray(x)
y = np.asarray(y)
r = np.sqrt(x**2 + y**2)  # radial position in the 2D Cartesian wedge
p_vals = np.asarray(p_vals)
T_vals = np.asarray(T_vals)

# Sort by radius
order = np.argsort(r)
r_sorted = r[order]
p_sorted = p_vals[order]
T_sorted = T_vals[order]
r_rel = r_sorted / Ro


def radial_bin_mean(r_rel_vals, field_vals, nbins=50):
    """Average a 2D field onto radial bins for a smooth 1D profile.

    The plot collapses a 2D wedge onto radius, so several nodes at different
    angles share one radius. Connecting them sorted by radius produces a
    sawtooth at the steep void front (angular spread reaches ~0.5 there). The
    bin mean removes that angular projection noise; the radial mean itself is
    monotone, so the curve is faithful, not just cosmetically smoothed.
    """
    edges = np.linspace(0.0, r_rel_vals.max(), nbins + 1)
    idx = np.digitize(r_rel_vals, edges)
    centres, means = [], []
    for b in range(1, nbins + 1):
        sel = idx == b
        if not np.any(sel):
            continue
        centres.append(0.5 * (edges[b - 1] + edges[b]))
        means.append(field_vals[sel].mean())
    return np.asarray(centres), np.asarray(means)


r_bin, p_bin = radial_bin_mean(r_rel, p_sorted)
_, T_bin = radial_bin_mean(r_rel, T_sorted)

# --. Derived quantities --..
VOID_THRESHOLD = 0.90          # p above this counts as central void
COLUMNAR_THRESHOLD = 0.02      # p below this counts as restructured (Barani Sec. 5.2)
RIM_FRACTION = 0.93            # outer band (r/Ro > this) used for the rim porosity

p_center = float(p_sorted[0])
T_center = float(T_sorted[0])
T_rim = float(T_sorted[-1])

# Central void radius: outermost radius still at void porosity.
void_nodes = r_rel[p_sorted >= VOID_THRESHOLD]
void_radius_rel = float(void_nodes.max()) if void_nodes.size > 0 else 0.0

# Rim porosity (genuine prediction, no imposed inflow).
rim_mask = r_rel > RIM_FRACTION
p_rim = float(p_sorted[rim_mask].mean()) if np.any(rim_mask) else float(p_sorted[-1])

# Restructured zone: does a low-porosity columnar band exist between void and rim?
p_min = float(p_sorted.min())
columnar_nodes = r_rel[(p_sorted < COLUMNAR_THRESHOLD) & (r_rel > void_radius_rel)]
columnar_outer_rel = float(columnar_nodes.max()) if columnar_nodes.size > 0 else void_radius_rel

print(f"[INFO] Centre temperature   : {T_center:.2f} K")
print(f"[INFO] Outer rim temperature: {T_rim:.2f} K")
print(f"[INFO] Centre porosity      : {p_center:.4f} (expected ~ 1.0)")
print(f"[INFO] Minimum porosity     : {p_min:.4f} (restructured zone)")
print(f"[INFO] Central void radius  : {void_radius_rel:.4f} r/Ro (Barani ~ 0.20)")
print(f"[INFO] Columnar zone edge   : {columnar_outer_rel:.4f} r/Ro")
print(f"[INFO] Rim porosity         : {p_rim:.4f} (expected ~ 0.15)")

# --. Paper-anchored absolute checks --..
VOID_RADIUS_REF = 0.20  # Barani Fig. 4: gradient/velocity null inside ~0.2 r/Ro
errors = {
    "center_void_porosity": {
        "numerical": p_center,
        "reference": 1.0,
        "abs_error": float(abs(p_center - 1.0)),
        "rel_error": float(abs(p_center - 1.0) / 1.0),
    },
    "void_radius_relative": {
        "numerical": void_radius_rel,
        "reference": VOID_RADIUS_REF,
        "abs_error": float(abs(void_radius_rel - VOID_RADIUS_REF)),
        "rel_error": float(abs(void_radius_rel - VOID_RADIUS_REF) / VOID_RADIUS_REF),
    },
    "rim_fabricated_porosity": {
        "numerical": p_rim,
        "reference": 0.15,
        "abs_error": float(abs(p_rim - 0.15)),
        "rel_error": float(abs(p_rim - 0.15) / 0.15),
    },
}

# Centre temperature: tracked for regression reproducibility only (no paper value).
temperature_info = {
    "center_temperature": {
        "numerical": T_center,
        "reference": T_center,
        "abs_error": 0.0,
        "rel_error": 0.0,
    }
}

# --. Plot 1: porosity radial profile --..
try:
    plt.figure(figsize=(7, 5))
    plt.scatter(r_rel, p_sorted, s=6, color="r", alpha=0.18, label="nodes (all angles)")
    plt.plot(r_bin, p_bin, "r-", lw=2.5, label="Z3ST porosity (radial mean)")
    plt.axhline(0.15, color="gray", ls="--", label="initial porosity (0.15)")
    plt.axvline(VOID_RADIUS_REF, color="k", ls=":", alpha=0.7, label="Barani void radius ~0.2")
    plt.xlabel("Relative radius r / Ro (-)")
    plt.ylabel("Porosity (-)")
    plt.title("Radial porosity profile at t = 10,000 s")
    plt.grid(True, ls=":", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "porosity_radial_profile.png"), dpi=150)
    print("[INFO] porosity_radial_profile.png saved")
except Exception as e:
    print(f"[WARNING] porosity plot skipped: {e}")

# --. Plot 2: temperature radial profile --..
try:
    plt.figure(figsize=(7, 5))
    plt.scatter(r_rel, T_sorted, s=6, color="b", alpha=0.18, label="nodes (all angles)")
    plt.plot(r_bin, T_bin, "b-", lw=2.5, label="Z3ST temperature (radial mean)")
    plt.xlabel("Relative radius r / Ro (-)")
    plt.ylabel("Temperature (K)")
    plt.title("Radial temperature profile at t = 10,000 s")
    plt.grid(True, ls=":", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "temperature_radial_profile.png"), dpi=150)
    print("[INFO] temperature_radial_profile.png saved")
except Exception as e:
    print(f"[WARNING] temperature plot skipped: {e}")

# --. Verdicts --..
# Tolerance is intentionally loose: the central-void radius depends on the pore
# velocity calibration (Sens constants), which is not yet pinned to a published
# table — see porosity_migration_model.py. This is a coarse engineering check.
all_checks = {**errors, **temperature_info}
# center_temperature carries rel_error 0 (reference == numerical), so it never
# drives the absolute pass/fail; it is written to the gold purely so the
# regression check tracks the centre temperature between runs.
pass_fail_check(all_checks, 0.30, OUT_JSON, CASE_DIR)
regression_check(all_checks, case_dir=CASE_DIR)

print("\n[INFO] Porosity migration non-regression completed.\n")
