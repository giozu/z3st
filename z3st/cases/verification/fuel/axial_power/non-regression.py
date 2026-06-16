#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: verification/fuel/axial_power

Thermo-mechanical verification of the axial power source bus
Five closed-form checks, independent of each other:

  source bus
  1. accumulation magnitude : bu_mean = q_avg * t_total / (rho * HM * 8.64e10),
     q_avg = lhr/area (the mean-1 normalisation preserves the average rating).
  2. axial peaking factor   : bu_peak/bu_mean = 1 / [(2 L'/pi L) sin(pi L/2 L')].
  3. end/peak shape         : bu(0)/bu(z_mid) = cos(pi L / 2 L').
  thermal
  4. centreline rise        : dT(z_mid) = q'''_peak R^2/(4k), q'''_peak = q_avg/F_mean.
  mechanical
  5. free elongation        : u_z(top) = alpha * LHR * L / (8 pi k)  (T_ref = T_cool).
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

from z3st.utils.utils_extract_vtu import extract_field
from z3st.utils.utils_verification import pass_fail_check, regression_check

CASE_DIR = os.path.dirname(__file__)
OUT = os.path.join(CASE_DIR, "output")
OUT_JSON = os.path.join(OUT, "non-regression.json")
TOLERANCE = 3e-2

_single = os.path.join(OUT, "fields.vtu")
_steps = sorted(glob.glob(os.path.join(OUT, "fields_*.vtu")))
VTU_FILE = _steps[-1] if _steps else _single

# --. geometry + material + history from the case YAML files --..
geom = yaml.safe_load(open(os.path.join(CASE_DIR, "geometry.yaml")))
Ro, L = float(geom["Ro"]), float(geom["Lz"])
area = np.pi * Ro**2

inp = yaml.safe_load(open(os.path.join(CASE_DIR, "input.yaml")))
lhr = float(inp["lhr"][-1])                    # constant linear heat rate (W/m)
t_total = float(max(inp["time"]))              # total irradiation time (s)

mat = yaml.safe_load(open(os.path.join(CASE_DIR, next(iter(inp["materials"].values())))))
rho = float(mat["rho"])
hm = float(mat.get("heavy_metal_fraction", 0.8815))
Lp = float(mat["axial_extrapolated_length"])
k = float(mat["k"]); alpha = float(mat["alpha"]); T_cool = float(mat["T_ref"])

# --. analytical references --..
SECONDS_PER_MWD = 8.64e10                       # 86400 s/day * 1e6 W/MW
q_avg = lhr / area
x = np.pi * L / (2.0 * Lp)
F_MEAN = np.sin(x) / x                          # continuous mean of cos profile
peaking = 1.0 / F_MEAN

BU_MEAN_REF = q_avg * t_total / (rho * hm * SECONDS_PER_MWD)   # 1
PEAK_RATIO_REF = peaking                                       # 2
END_RATIO_REF = float(np.cos(x))                              # 3
DT_MID_REF = (q_avg * peaking) * Ro**2 / (4.0 * k)            # 4
UZ_TOP_REF = alpha * lhr * L / (8.0 * np.pi * k)             # 5
P_REF = lhr * L

# --. numerical results --..
xr, yz, _, bu = extract_field(VTU_FILE, field_name="Burnup")
z = np.asarray(yz)                              # axisymmetric (r-z) mesh: y = z
bu = np.asarray(bu)
bu_mean = float(np.mean(bu))
bu_peak = float(bu[np.argmin(np.abs(z - 0.5 * L))])   # node at z_mid
bu_end = float(bu[np.argmin(z)])                       # node at z = 0
peak_ratio = bu_peak / max(bu_mean, 1e-12)
end_ratio = bu_end / max(bu_peak, 1e-12)

# temperature + displacement (vector) — read directly for the r-z components
g = pv.read(VTU_FILE)
r_pt, z_pt = g.points[:, 0], g.points[:, 1]
T = np.asarray(g.point_data["Temperature"]).ravel()
uz = np.asarray(g.point_data["Displacement"])[:, 1]   # axial component (r-z mesh)

axis = r_pt < 1e-6
za, dTa = z_pt[axis], T[axis] - T_cool
o = np.argsort(za); za, dTa = za[o], dTa[o]
dT_mid = float(dTa[np.argmin(np.abs(za - 0.5 * L))])
uz_top = float(uz[z_pt > (L - 1e-6)].mean())

print(f"[INFO] burnup mean     : {bu_mean:.6e}  (ref {BU_MEAN_REF:.6e} MWd/kgU)")
print(f"[INFO] axial peaking   : {peak_ratio:.4f}      (ref {PEAK_RATIO_REF:.4f})")
print(f"[INFO] end/peak ratio  : {end_ratio:.4f}      (ref {END_RATIO_REF:.4f})")
print(f"[INFO] dT centre @z_mid: {dT_mid:.2f} K   (ref {DT_MID_REF:.2f})")
print(f"[INFO] elongation u_top: {uz_top*1e6:.2f} um  (ref {UZ_TOP_REF*1e6:.2f})")


def _entry(num, ref):
    return {"numerical": float(num), "reference": float(ref),
            "abs_error": float(abs(num - ref)), "rel_error": float(abs(num - ref) / abs(ref))}


errors = {
    "burnup_mean_closed_form": _entry(bu_mean, BU_MEAN_REF),
    "axial_peaking_factor": _entry(peak_ratio, PEAK_RATIO_REF),
    "end_peak_ratio": _entry(end_ratio, END_RATIO_REF),
    "dT_centre_mid": _entry(dT_mid, DT_MID_REF),
    "elongation_u_top": _entry(uz_top, UZ_TOP_REF),
}

# integrated power from the set_power diagnostic (separable shaping -> LHR*L)
LOG = os.path.join(CASE_DIR, "log_z3st.md")
if os.path.exists(LOG):
    hits = re.findall(r"Integrated fissile power in \S+:\s*([0-9.eE+\-]+)", open(LOG).read())
    if hits:
        errors["total_power"] = _entry(float(hits[-1]), P_REF)

# --. figures --..
z_line = np.linspace(0.0, L, 200)
f_line = np.cos(np.pi * (z_line - 0.5 * L) / Lp) / F_MEAN
plt.figure(figsize=(7, 5))
plt.plot(BU_MEAN_REF * f_line, z_line * 1e2, "k-", lw=2.5, alpha=0.7,
         label="analytical  bu_mean·cos(π(z−z_mid)/L′)/⟨f⟩")
plt.scatter(bu, z * 1e2, s=10, facecolors="none", edgecolors="r", label="Z3ST nodal burnup")
plt.axvline(BU_MEAN_REF, color="0.5", ls=":", lw=1.2, label=f"mean = {BU_MEAN_REF:.1f} MWd/kgU")
plt.xlabel("burnup (MWd/kgU)"); plt.ylabel("axial position z (cm)")
plt.title(f"Axial burnup: chopped cosine (L={L:.2f} m, L′={Lp:.2f} m, peaking={peaking:.3f})")
plt.grid(True, ls=":", alpha=0.6); plt.legend(); plt.tight_layout()
plt.savefig(os.path.join(OUT, "burnup_axial_profile.png"), dpi=150)
print("[INFO] burnup_axial_profile.png saved")

plt.figure(figsize=(7, 5))
plt.plot(DT_MID_REF * np.cos(np.pi * (z_line - 0.5 * L) / Lp), z_line * 1e2, "k-", lw=2.0,
         alpha=0.7, label=r"analytic $q'''(z)R^2/4k$")
plt.scatter(dTa, za * 1e2, s=12, facecolors="none", edgecolors="C0", label="Z3ST centreline rise")
plt.xlabel(r"$T_\mathrm{centre}-T_\mathrm{cool}$ (K)"); plt.ylabel("axial position z (cm)")
plt.title("Axial centreline temperature rise"); plt.grid(True, ls=":", alpha=0.6)
plt.legend(); plt.tight_layout()
plt.savefig(os.path.join(OUT, "temperature_axial_profile.png"), dpi=150)
print("[INFO] temperature_axial_profile.png saved")

pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)
print("\n[INFO] non-regression completed.\n")
