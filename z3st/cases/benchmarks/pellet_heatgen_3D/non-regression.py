#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: pellet_heatgen_3D

Solid UO2 pellet (diameter 1 cm) in 3D at ~25 kW/m, rim at 600 K. Half model:
5 mm meshed with a symmetry plane (u_z = 0) on the bottom face. Constant k keeps
the analytic radial temperature profile valid; the script also checks for
spurious (hourglass-like) modes:

  - axisymmetry: azimuthal scatter of T and von Mises at fixed radius;
  - orthocylindricity: in-plane shear sigma_r-theta vs peak hoop;
  - bulk response: peak temperature, displacement, von Mises.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import yaml

from z3st.utils.utils_extract_xdmf import *
from z3st.utils.utils_verification import *

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
XDMF_FILE = os.path.join(CASE_DIR, "output", "fields.xdmf")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")

with open(os.path.join(CASE_DIR, "geometry.yaml")) as f:
    geom = yaml.safe_load(f)
Ri, Ro, Lz = float(geom["Ri"]), float(geom["Ro"]), float(geom["Lz"])  # m

with open(os.path.join(CASE_DIR, "input.yaml")) as f:
    inp = yaml.safe_load(f)
mat_path = os.path.join(CASE_DIR, next(iter(inp["materials"].values())))
with open(mat_path) as f:
    mat = yaml.safe_load(f)
LHR = float(inp["lhr"][0])  # W/m

k = float(mat["k"])  # W/m.K (constant)
E, nu, alpha = float(mat["E"]), float(mat["nu"]), float(mat["alpha"])
To = 600.0  # K   rim temperature (boundary_conditions.yaml)

# Extract on the symmetry plane (z = 0): mid-pellet, closest to the
# infinite-cylinder plane-strain state.
n_z = 12  # extrusion layers (mesh.geo)
z_target = 0.0
z_tol = Lz / n_z

# Loose stress tolerance: a 1 cm pellet is too short for plane strain, so the
# infinite-cylinder analytic is only ~30 % accurate. The tight checks are the
# analytic temperature and the regression vs gold; stress is a diagnostic only.
TOLERANCE = 3.5e-1


# --.. ..- .-.. .-.. --- analytic temperature (constant k) --.. ..- .-.. .-.. ---
def analytic_T(r):
    r = np.asarray(r)
    return LHR / (4 * np.pi * k) * (1 - r**2 / Ro**2) + To


def analytical_thermal_stress(x, T):
    """Plane-strain pipe-equation thermal stress (reference shape)."""
    sigma_star = alpha * E * (max(T) - min(T)) / (4 * (1.0 - nu))
    sigma_tt = -sigma_star * (1.0 - 3.0 * (x / Ro) ** 2)
    sigma_rr = -sigma_star * (1.0 - (x / Ro) ** 2)
    sigma_zz = sigma_tt + sigma_rr
    return sigma_rr, sigma_tt, sigma_zz


def to_cylindrical(S, xc, yc):
    """Rotate flattened 3x3 Cartesian stress (N,9) into (r, theta, z)."""
    r = np.sqrt(xc**2 + yc**2)
    r_safe = np.where(r > 0, r, 1.0)
    c, s = xc / r_safe, yc / r_safe
    sxx, sxy, syy, szz = S[:, 0], S[:, 1], S[:, 4], S[:, 8]
    s_rr = c * c * sxx + s * s * syy + 2.0 * c * s * sxy
    s_tt = s * s * sxx + c * c * syy - 2.0 * c * s * sxy
    s_rt = -c * s * sxx + c * s * syy + (c * c - s * s) * sxy
    return s_rr, s_tt, szz, s_rt


def von_mises(S):
    sxx, syy, szz = S[:, 0], S[:, 4], S[:, 8]
    sxy, syz, sxz = S[:, 1], S[:, 5], S[:, 2]
    return np.sqrt(
        0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2)
        + 3.0 * (sxy**2 + syz**2 + sxz**2)
    )


def azimuthal_scatter(values, radii, scale, nbins=20):
    """Max over radial bins of the azimuthal std, normalised by ``scale``."""
    edges = np.linspace(0.0, Ro, nbins + 1)
    which = np.clip(np.digitize(radii, edges) - 1, 0, nbins - 1)
    worst = 0.0
    for b in range(nbins):
        v = values[which == b]
        if v.size > 2:
            worst = max(worst, float(np.std(v)))
    return worst / max(scale, 1e-30)


# --.. ..- .-.. .-.. --- read fields --.. ..- .-.. .-.. ---
list_fields_xdmf(XDMF_FILE)
print(f"[INFO] Mid-height extraction plane: z = {z_target:.4e} m (tol {z_tol:.2e})")

# Temperature (nodal)
xT, yT, zT, T_all = extract_field_xdmf(XDMF_FILE, field_name="Temperature", step_index=-1)
rT_full = np.sqrt(xT**2 + yT**2)
mT = np.abs(zT - z_target) < z_tol
sortT = np.argsort(rT_full[mT])
r_T = rT_full[mT][sortT]
T = T_all[mT][sortT]
T_ref = analytic_T(r_T)

# Displacement (nodal) -- peak magnitude over the whole pellet
xU, yU, zU, U_all = extract_field_xdmf(XDMF_FILE, field_name="Displacement", step_index=-1)
u_mag = np.sqrt((U_all**2).sum(axis=1))
u_max = float(np.max(u_mag))

# Stress (cell data, Cartesian)
xS, yS, zS, S_all = extract_field_xdmf(XDMF_FILE, field_name="Stress", step_index=-1)
mS = np.abs(zS - z_target) < z_tol
xc, yc = xS[mS], yS[mS]
S_mid = S_all[mS]
r_s = np.sqrt(xc**2 + yc**2)
s_rr, s_tt, s_zz, s_rt = to_cylindrical(S_mid, xc, yc)
vm = von_mises(S_mid)
vm_all = von_mises(S_all)
vm_max = float(np.max(vm_all))

sortS = np.argsort(r_s)
r_s_s = r_s[sortS]
s_rr_s, s_tt_s, s_zz_s = s_rr[sortS], s_tt[sortS], s_zz[sortS]
s_rr_ana, s_tt_ana, s_zz_ana = analytical_thermal_stress(r_s_s, analytic_T(r_s_s))

# --.. ..- .-.. .-.. --- spurious-mode diagnostics --.. ..- .-.. .-.. ---
dT_scale = max(np.max(T_ref) - np.min(T_ref), 1e-12)
azim_T = azimuthal_scatter(T, r_T, dT_scale)
azim_vm = azimuthal_scatter(vm, r_s, max(np.mean(vm), 1e-30))
ortho_shear = float(np.max(np.abs(s_rt)) / max(np.max(np.abs(s_tt)), 1e-30))

print(f"[INFO] Peak temperature      : {np.max(T) - 273.15:.1f} C ({np.max(T):.1f} K)")
print(f"[INFO] Peak displacement     : {u_max * 1e6:.2f} um")
print(f"[INFO] Peak von Mises stress : {vm_max * 1e-6:.1f} MPa")
print(f"[INFO] Azimuthal T scatter   : {azim_T:.3e}")
print(f"[INFO] Azimuthal vM scatter  : {azim_vm:.3e}  (hourglass / checkerboard indicator)")
print(f"[INFO] Orthocylindricity shear max|s_rt|/max|s_tt|: {ortho_shear:.3e}")

# --.. ..- .-.. .-.. --- plot --.. ..- .-.. .-.. ---
Pa_to_MPa = 1e-6
def radial_profile(values, radii, nbins=24):
    """Azimuthal mean +/- std of ``values`` over ``nbins`` radial bins."""
    edges = np.linspace(0.0, Ro, nbins + 1)
    rc = 0.5 * (edges[:-1] + edges[1:])
    which = np.clip(np.digitize(radii, edges) - 1, 0, nbins - 1)
    mean = np.full(nbins, np.nan)
    std = np.zeros(nbins)
    for b in range(nbins):
        v = values[which == b]
        if v.size:
            mean[b] = v.mean()
            std[b] = v.std()
    ok = ~np.isnan(mean)
    return rc[ok], mean[ok], std[ok]


plt.figure(figsize=(10, 7))
ax1 = plt.gca()
for vals, ana, col, lab in [
    (s_rr, s_rr_ana, "r", r"\sigma_{rr}"),
    (s_tt, s_tt_ana, "g", r"\sigma_{\theta\theta}"),
    (s_zz, s_zz_ana, "b", r"\sigma_{zz}"),
]:
    rc, mu, sd = radial_profile(vals, r_s)
    ax1.fill_between(rc, (mu - sd) * Pa_to_MPa, (mu + sd) * Pa_to_MPa, color=col, alpha=0.15)
    ax1.plot(rc, mu * Pa_to_MPa, col + "o", ms=4, label=rf"Num. ${lab}$ (azim. mean)")
ax1.plot(r_s_s, s_rr_ana * Pa_to_MPa, "r-", lw=1.5, label=r"Ana. $\sigma_{rr}$")
ax1.plot(r_s_s, s_tt_ana * Pa_to_MPa, "g-", lw=1.5, label=r"Ana. $\sigma_{\theta\theta}$")
ax1.plot(r_s_s, s_zz_ana * Pa_to_MPa, "b-", lw=1.5, label=r"Ana. $\sigma_{zz}$")
ax1.set_xlabel("Radius (m)")
ax1.set_ylabel("Stress (MPa)")
ax1.grid(True, ls="--", alpha=0.7)
ax2 = ax1.twinx()
rcT, muT, _ = radial_profile(T, r_T)
ax2.plot(rcT, muT - 273.15, "ks", ms=3, alpha=0.6, label="Num. T (azim. mean)")
ax2.plot(r_T, T_ref - 273.15, "k--", lw=1.0, alpha=0.8, label="Ana. T")
ax2.set_ylabel("Temperature (C)")
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines + lines2, labels + labels2, loc="best", frameon=True)
plt.title(
    rf"UO2 pellet, LHR = {LHR/1e3:.0f} kW/m: $T_{{max}}$ = {np.max(T)-273.15:.0f} C, "
    rf"$\sigma_{{vM,max}}$ = {vm_max*1e-6:.0f} MPa, azim. vM scatter = {azim_vm:.1e}",
    pad=15,
)
plt.tight_layout()
plot_path = os.path.join(CASE_DIR, "output", "stress_comparison.png")
plt.savefig(plot_path, dpi=200)
print(f"[INFO] Plot saved in: {plot_path}")


# --.. ..- .-.. .-.. --- metrics --.. ..- .-.. .-.. ---
def rel_l2(num, ref):
    return float(np.sqrt(np.mean((num - ref) ** 2)) / np.sqrt(np.mean(ref**2)))


L2_T = float(np.sqrt(np.mean((T - T_ref) ** 2)))
RelL2_T = float(L2_T / np.mean(np.abs(T_ref)))
Tmax_num, Tmax_ref = float(np.max(T)), float(np.max(T_ref))

errors = {
    "L2_error_T": {"numerical": L2_T, "reference": 0.0, "abs_error": L2_T, "rel_error": RelL2_T},
    "T_max": {
        "numerical": Tmax_num, "reference": Tmax_ref,
        "abs_error": abs(Tmax_num - Tmax_ref), "rel_error": abs(Tmax_num - Tmax_ref) / Tmax_ref,
    },
    "u_max": {"numerical": u_max, "reference": u_max, "abs_error": 0.0, "rel_error": 0.0},
    "vonMises_max": {"numerical": vm_max, "reference": vm_max, "abs_error": 0.0, "rel_error": 0.0},
    "azimuthal_nonuniformity_T": {
        "numerical": azim_T, "reference": 0.0, "abs_error": azim_T, "rel_error": azim_T,
    },
    "azimuthal_nonuniformity_vonMises": {
        "numerical": azim_vm, "reference": 0.0, "abs_error": azim_vm, "rel_error": azim_vm,
    },
    "orthocylindricity_shear": {
        "numerical": ortho_shear, "reference": 0.0, "abs_error": ortho_shear, "rel_error": ortho_shear,
    },
    "L2_error_sigma_rr": {"numerical": rel_l2(s_rr_s, s_rr_ana), "reference": 0.0,
                          "abs_error": rel_l2(s_rr_s, s_rr_ana), "rel_error": rel_l2(s_rr_s, s_rr_ana)},
    "L2_error_sigma_tt": {"numerical": rel_l2(s_tt_s, s_tt_ana), "reference": 0.0,
                          "abs_error": rel_l2(s_tt_s, s_tt_ana), "rel_error": rel_l2(s_tt_s, s_tt_ana)},
    "L2_error_sigma_zz": {"numerical": rel_l2(s_zz_s, s_zz_ana), "reference": 0.0,
                          "abs_error": rel_l2(s_zz_s, s_zz_ana), "rel_error": rel_l2(s_zz_s, s_zz_ana)},
}

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] Non-regression completed.\n")
