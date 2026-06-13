#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: verification/mechanics/tall_rod_axial_2D

Thermo-mechanical verification of the axial power form factor on a slender rod (R << L): 
Because the rod is slender, conduction is locally 1-D radial at each elevation
  1. integrated power     P = LHR * L
     (the mean-1 normalisation of the axial factor preserves the rating).
  2. centreline rise      dT(z_mid) = q'''_peak * R^2 / (4 k),
     with q'''_peak = (LHR/A) / F_mean  the mid-height volumetric rating.
  3. axial peaking        peak/mean of the centreline rise = 1 / F_mean,
     F_mean = (2 L' / pi L) sin(pi L / 2 L')  the cosine's continuous mean.
  4. free elongation      u_z(top) = alpha * LHR * L / (8 pi k)

"""

import os
import re
import glob
import yaml
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pyvista as pv

from z3st.utils.utils_verification import pass_fail_check, regression_check

CASE = os.path.dirname(__file__)
OUT = os.path.join(CASE, "output")
OUT_JSON = os.path.join(OUT, "non-regression.json")
TOLERANCE = 3.0e-2

# --. case parameters from the YAML files --..
geom = yaml.safe_load(open(os.path.join(CASE, "geometry.yaml")))
Ro, L = float(geom["Ro"]), float(geom["Lz"])
area = np.pi * Ro**2

inp = yaml.safe_load(open(os.path.join(CASE, "input.yaml")))
lhr = float(inp["lhr"][-1])
mat = yaml.safe_load(open(os.path.join(CASE, next(iter(inp["materials"].values())))))
k = float(mat["k"]); alpha = float(mat["alpha"])
T_cool = float(mat["T_ref"])
Lp = float(mat["axial_extrapolated_length"])

# --. analytical references --..
x = np.pi * L / (2.0 * Lp)
F_mean = np.sin(x) / x                       # continuous mean of the cosine
peaking = 1.0 / F_mean
q_avg = lhr / area
q_peak = q_avg * peaking
P_REF = lhr * L                              # 1. integrated power (W)
DT_MID_REF = q_peak * Ro**2 / (4.0 * k)      # 2. centreline rise at z_mid (K)
PEAK_REF = peaking                           # 3. axial peaking of dT
UZ_TOP_REF = alpha * lhr * L / (8.0 * np.pi * k)   # 4. free elongation (m)

# --. numerical results --..
files = sorted(glob.glob(os.path.join(OUT, "fields_*.vtu")))
VTU = files[-1] if files else os.path.join(OUT, "fields.vtu")
m = pv.read(VTU)
r, z = m.points[:, 0], m.points[:, 1]
T = np.asarray(m.point_data["Temperature"]).ravel()
uz = np.asarray(m.point_data["Displacement"])[:, 1]   # axial component (r-z mesh)

# centreline (axis) temperature rise vs z
axis = r < 1e-6
za, dTa = z[axis], T[axis] - T_cool
o = np.argsort(za); za, dTa = za[o], dTa[o]
dT_mid = float(dTa[np.argmin(np.abs(za - 0.5 * L))])
peak_num = float(dTa.max() / max(dTa.mean(), 1e-12))

# mean axial displacement of the top face
top = z > (L - 1e-6)
uz_top = float(uz[top].mean())

# integrated power from the set_power diagnostic in the run log
P_num = np.nan
log = os.path.join(CASE, "log_z3st.md")
if os.path.exists(log):
    mlog = re.findall(r"Integrated fissile power in \w+:\s*([-\d.eE+]+)", open(log).read())
    if mlog:
        P_num = float(mlog[-1])

print(f"[INFO] integrated power : {P_num:.4e} W   (ref {P_REF:.4e})")
print(f"[INFO] dT centre @z_mid : {dT_mid:.3f} K   (ref {DT_MID_REF:.3f})")
print(f"[INFO] axial peaking    : {peak_num:.4f}     (ref {PEAK_REF:.4f})")
print(f"[INFO] elongation u_top : {uz_top*1e6:.3f} um  (ref {UZ_TOP_REF*1e6:.3f})")

# --. axial profile figure --..
fig, ax = plt.subplots(figsize=(6.5, 4.5))
ax.plot(za * 1e3, dTa, "o", ms=3, color="#C44E52", label="Z3ST centreline rise")
ax.plot(za * 1e3, DT_MID_REF * np.cos(np.pi * (za - 0.5 * L) / Lp), "k--", lw=1.4,
        label=r"analytic $q'''(z)R^2/4k$")
ax.set_xlabel("axial position z (mm)"); ax.set_ylabel(r"$T_\mathrm{centre}-T_\mathrm{cool}$ (K)")
ax.set_title("Axial centreline temperature rise"); ax.grid(alpha=0.3); ax.legend()
fig.tight_layout(); fig.savefig(os.path.join(OUT, "axial_profile.png"), dpi=150)
print("[INFO] axial_profile.png saved")


def _entry(num, ref):
    return {"numerical": num, "reference": ref,
            "abs_error": abs(num - ref), "rel_error": abs(num - ref) / abs(ref)}


errors = {
    "integrated_power": _entry(P_num, P_REF),
    "dT_centre_mid": _entry(dT_mid, DT_MID_REF),
    "axial_peaking": _entry(peak_num, PEAK_REF),
    "elongation_u_top": _entry(uz_top, UZ_TOP_REF),
}

pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE)
regression_check(errors, CASE)
print("\n[INFO] non-regression completed.\n")
