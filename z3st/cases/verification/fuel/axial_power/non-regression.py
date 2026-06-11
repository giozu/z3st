#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: verification/fuel/axial_power

Verifies the AXIAL power source bus (axial_profile -> set_power f(z)) on a tall
axisymmetric fuel column — the Todreas & Kazimi "1-D axial problem".

A constant linear heat rate is shaped by the chopped-cosine form factor

    f(z) = cos(pi (z - z_mid) / L'),

with L the active (meshed) height and L' the extrapolated length, normalised to
nodal mean 1 by set_power. Burnup accumulates as bu(z) = q(z)·t/(rho·HM·8.64e10),
so the final burnup field IS the normalised profile, and three closed-form
checks follow, independent of each other:

  1. accumulation magnitude — the nodal-mean burnup equals the flat closed form
         bu_mean = q_avg * t_total / (rho * HM * 8.64e10),
     q_avg = lhr / area: the mean-1 normalisation preserves the average rating.
  2. axial peaking factor — peak/mean burnup equals
         1 / [ (2 L' / pi L) sin(pi L / 2 L') ],
     the inverse of the continuous mean of the un-normalised cosine.
  3. end/peak shape — bu(z=0) / bu(z_mid) = cos(pi L / 2 L'), independent of
     the normalisation.

One figure: the axial burnup profile, Z3ST nodal vs analytical.
"""

import os
import re
import glob
import yaml
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_verification import *

CASE_DIR = os.path.dirname(__file__)
OUT = os.path.join(CASE_DIR, "output")
OUT_JSON = os.path.join(OUT, "non-regression.json")

# Burnup accumulates over the run -> read the *final* step.
_single = os.path.join(OUT, "fields.vtu")
_steps = sorted(glob.glob(os.path.join(OUT, "fields_*.vtu")))
VTU_FILE = _steps[-1] if _steps else _single

# --. geometry + material + history from the case YAML files --..
with open(os.path.join(CASE_DIR, "geometry.yaml")) as f:
    geom = yaml.safe_load(f)
Ro = float(geom["Ro"])
L = float(geom["Lz"])                          # active fuel height
area = np.pi * Ro**2

with open(os.path.join(CASE_DIR, "input.yaml")) as f:
    inp = yaml.safe_load(f)
lhr = float(inp["lhr"][-1])                    # constant linear heat rate (W/m)
t_total = float(max(inp["time"]))              # total irradiation time (s)

with open(os.path.join(CASE_DIR, next(iter(inp["materials"].values())))) as f:
    mat = yaml.safe_load(f)
rho = float(mat["rho"])
hm = float(mat.get("heavy_metal_fraction", 0.8815))
L_prime = float(mat["axial_extrapolated_length"])

# --. analytical references --..
SECONDS_PER_MWD = 8.64e10                      # 86400 s/day * 1e6 W/MW
q_avg = lhr / area
BU_MEAN_REF = q_avg * t_total / (rho * hm * SECONDS_PER_MWD)

x = np.pi * L / (2.0 * L_prime)                # pi L / 2 L'
F_MEAN = np.sin(x) / x                         # continuous mean of cos profile
PEAK_RATIO_REF = 1.0 / F_MEAN                  # axial peaking factor
END_RATIO_REF = float(np.cos(x))               # f(end) / f(mid)
TOLERANCE = 1e-2

# --. numerical results --..
xr, yz, _, bu = extract_field(VTU_FILE, field_name="Burnup")
z = np.asarray(yz)                             # axisymmetric (r-z) mesh: y = z
bu = np.asarray(bu)

bu_mean = float(np.mean(bu))
bu_peak = float(bu[np.argmin(np.abs(z - 0.5 * L))])   # node at z_mid
bu_end = float(bu[np.argmin(z)])                       # node at z = 0
peak_ratio = bu_peak / max(bu_mean, 1e-12)
end_ratio = bu_end / max(bu_peak, 1e-12)

print(f"[INFO] q_avg = lhr/area = {q_avg:.4e} W/m^3 over t = {t_total:.3e} s")
print(f"[INFO] nodal-mean burnup : numerical = {bu_mean:.6e}, "
      f"analytical = {BU_MEAN_REF:.6e} MWd/kgU")
print(f"[INFO] peak/mean ratio   : numerical = {peak_ratio:.4f}, "
      f"analytical = {PEAK_RATIO_REF:.4f}")
print(f"[INFO] end/peak ratio    : numerical = {end_ratio:.4f}, "
      f"analytical cos(piL/2L') = {END_RATIO_REF:.4f}")

errors = {
    "burnup_mean_closed_form": {
        "numerical": bu_mean,
        "reference": BU_MEAN_REF,
        "abs_error": float(abs(bu_mean - BU_MEAN_REF)),
        "rel_error": float(abs(bu_mean - BU_MEAN_REF) / BU_MEAN_REF),
    },
    "axial_peaking_factor": {
        "numerical": peak_ratio,
        "reference": PEAK_RATIO_REF,
        "abs_error": float(abs(peak_ratio - PEAK_RATIO_REF)),
        "rel_error": float(abs(peak_ratio - PEAK_RATIO_REF) / PEAK_RATIO_REF),
    },
    "end_peak_ratio": {
        "numerical": end_ratio,
        "reference": END_RATIO_REF,
        "abs_error": float(abs(end_ratio - END_RATIO_REF)),
        "rel_error": float(abs(end_ratio - END_RATIO_REF) / END_RATIO_REF),
    },
}

# --. integrated power (parsed from the solver log) --..
# set_power prints the exact FE integral of the fissile source per step; the
# axial shaping is separable from the 2πr weight, so the integral must equal
# LHR·L (the normalisation preserves the rating).
LOG = os.path.join(CASE_DIR, "log_z3st.md")
P_REF = lhr * L
if os.path.exists(LOG):
    with open(LOG) as f:
        hits = re.findall(r"Integrated fissile power in \S+:\s*([0-9.eE+\-]+)", f.read())
    if hits:
        P_int = float(hits[-1])
        print(f"[INFO] integrated power  : numerical = {P_int:.6e} W, "
              f"analytical LHR·L = {P_REF:.6e} W")
        errors["total_power"] = {
            "numerical": P_int,
            "reference": P_REF,
            "abs_error": float(abs(P_int - P_REF)),
            "rel_error": float(abs(P_int - P_REF) / P_REF),
        }
    else:
        print("[WARNING] no 'Integrated fissile power' line in log — power check skipped.")
else:
    print("[WARNING] log_z3st.md not found — power check skipped.")

# --. figure: axial burnup profile (Z3ST nodal vs analytical) --..
try:
    z_line = np.linspace(0.0, L, 200)
    f_line = np.cos(np.pi * (z_line - 0.5 * L) / L_prime) / F_MEAN
    bu_line = BU_MEAN_REF * f_line
    order = np.argsort(z)
    plt.figure(figsize=(7, 5))
    plt.plot(bu_line, z_line * 1e2, "k-", lw=2.5, alpha=0.7,
             label="Analytical  bu_mean · cos(π(z−z_mid)/L′)/⟨f⟩")
    plt.scatter(bu, z * 1e2, s=10, facecolors="none", edgecolors="r",
                label="Z3ST nodal burnup")
    plt.axvline(BU_MEAN_REF, color="0.5", ls=":", lw=1.2,
                label=f"mean = {BU_MEAN_REF:.1f} MWd/kgU")
    plt.xlabel("burnup (MWd/kgU)")
    plt.ylabel("axial position z (cm)")
    plt.title("Axial burnup distribution: chopped-cosine form factor\n"
              f"L = {L:.2f} m, L′ = {L_prime:.2f} m, peaking = {PEAK_RATIO_REF:.3f}")
    plt.grid(True, ls=":", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "burnup_axial_profile.png"), dpi=150)
    print("[INFO] burnup_axial_profile.png saved")
except Exception as e:
    print(f"[WARNING] axial-profile plot skipped: {type(e).__name__}: {e}")

pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] non-regression completed.\n")
