#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: verification/mechanics/thermoelastic_cylinder_invariants_2D

Open-ended thick cylinder under internal pressure plus a radial temperature
gradient. Two checks:

1. the total stress against Lame (open end, sigma_zz = 0) superposed on the
   free-end thermal stresses: this is the one that tests the FE solution;
2. the split sigma = s + p I (deviatoric + hydrostatic): von Mises and Tresca
   of s equal those of sigma, and both vanish on p I. These are algebraic
   identities (both criteria depend only on differences of principal
   stresses), so they hold to round-off whatever the FE does: they test the
   tensor layout and the writer's VonMises / Tresca fields, not the physics.
"""

import os

import matplotlib.pyplot as plt
import numpy as np

from z3st.utils.non_regression import case_paths, error_metric, finish, load_case, metric
from z3st.utils.utils_extract_vtu import extract_field, list_fields

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR, VTU_FILE, OUT_JSON = case_paths(__file__)

# Geometry and material
geom, inp, mat = load_case(CASE_DIR)
Ri, Ro, Lz = float(geom["Ri"]), float(geom["Ro"]), float(geom["Lz"])  # m   inner/outer radius, height

E, nu, alpha = float(mat["E"]), float(mat["nu"]), float(mat["alpha"])  # Pa, -, 1/K

Pi, Po = 50.0e6, 0.0  # Pa         internal and external pressure
Ti, To = 500.0, 400.0  # K          inner and outer surface temperature
z_target, z_tol = Lz / 2, Lz / 10  # m          z-plane for data extraction

TOLERANCE = 3.0e-2  # -          tolerance for non-regression
ROUNDOFF = 1.0e-12  # -          tolerance for the algebraic identities


# --.. ..- .-.. .-.. --- analytic functions  --.. ..- .-.. .-.. ---
def analytic_T(r):
    """Logarithmic temperature profile with imposed T_inner and T_outer."""
    return Ti + (To - Ti) * np.log(r / Ri) / np.log(Ro / Ri)


def analytical_thermal_stress(r):
    """Free-end thermal stresses in a cylinder, integral form, any T(r)."""
    r_star = np.linspace(Ri, Ro, 20000)
    T_star = analytic_T(r_star)
    I_global = np.trapezoid(T_star * r_star, r_star)
    c = alpha * E / (1 - nu)

    sigma_r, sigma_t, sigma_z = np.zeros_like(r), np.zeros_like(r), np.zeros_like(r)
    for i, ri in enumerate(r):
        mask = r_star <= ri
        I_local = np.trapezoid(T_star[mask] * r_star[mask], r_star[mask])
        sigma_r[i] = c / ri**2 * (((ri**2 - Ri**2) / (Ro**2 - Ri**2)) * I_global - I_local)
        sigma_t[i] = c * (((1 + (Ri / ri) ** 2) / (Ro**2 - Ri**2)) * I_global + I_local / ri**2 - analytic_T(ri))
        sigma_z[i] = c * ((2 / (Ro**2 - Ri**2)) * I_global - analytic_T(ri))
    return sigma_r, sigma_t, sigma_z


def analytical_lame_open(r):
    """Lame stresses, open-ended cylinder (no axial load: sigma_zz = 0)."""
    A = (Pi * Ri**2 - Po * Ro**2) / (Ro**2 - Ri**2)
    B = (Ri**2 * Ro**2 * (Pi - Po)) / (Ro**2 - Ri**2)
    return A - B / r**2, A + B / r**2, np.zeros_like(r)


def von_mises(sig):
    """From the full (..., 3, 3) tensor, independent of the writer's formula."""
    s = sig - np.trace(sig, axis1=-2, axis2=-1)[..., None, None] / 3.0 * np.eye(3)
    return np.sqrt(1.5 * np.einsum("...ij,...ij", s, s))


def tresca(sig):
    """sigma_1 - sigma_3 from the principal stresses."""
    ev = np.linalg.eigvalsh(sig)
    return ev[..., -1] - ev[..., 0]


# --.. ..- .-.. .-.. --- checks --.. ..- .-.. .-.. ---
list_fields(VTU_FILE)

# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
# Whole domain: the identities must hold in every cell, end effects included
x_S, z_S, _, S_all = extract_field(VTU_FILE, field_name="Stress (cells)")
_, _, _, vm_writer = extract_field(VTU_FILE, field_name="VonMises (cells)")
_, _, _, tr_writer = extract_field(VTU_FILE, field_name="Tresca (cells)")

sig = S_all.reshape(-1, 3, 3)  # axisymmetric layout (r, theta, z)
p = np.trace(sig, axis1=-2, axis2=-1) / 3.0
sig_hyd = p[:, None, None] * np.eye(3)
sig_dev = sig - sig_hyd

vm_tot, vm_dev, vm_hyd = von_mises(sig), von_mises(sig_dev), von_mises(sig_hyd)
tr_tot, tr_dev, tr_hyd = tresca(sig), tresca(sig_dev), tresca(sig_hyd)
scale = np.max(vm_tot)

print(f"[INFO] max von Mises = {scale * 1e-6:.3f} MPa, max Tresca = {np.max(tr_tot) * 1e-6:.3f} MPa")
print(f"[INFO] max |p|        = {np.max(np.abs(p)) * 1e-6:.3f} MPa")

identities = {
    "vM_dev_minus_vM_tot": np.max(np.abs(vm_dev - vm_tot)) / scale,
    "vM_hyd": np.max(vm_hyd) / scale,
    "Tresca_dev_minus_Tresca_tot": np.max(np.abs(tr_dev - tr_tot)) / scale,
    "Tresca_hyd": np.max(tr_hyd) / scale,
    "writer_vM_minus_vM_tot": np.max(np.abs(vm_writer - vm_tot)) / scale,
    "writer_Tresca_minus_Tresca_tot": np.max(np.abs(tr_writer - tr_tot)) / scale,
}
for k, v in identities.items():
    print(f"[INFO] {k:32s} = {v:.3e}")
    assert v < ROUNDOFF, f"{k} = {v:.3e} exceeds round-off ({ROUNDOFF:.0e})"

# Bounds between the criteria: sqrt(3)/2 Tresca <= vM <= Tresca
assert np.all(vm_tot <= tr_tot * (1 + ROUNDOFF)), "von Mises exceeds Tresca"
assert np.all(vm_tot >= np.sqrt(3) / 2 * tr_tot * (1 - ROUNDOFF)), "von Mises below sqrt(3)/2 Tresca"

# Mid-plane profile against Lame + thermal
mid = np.abs(z_S - z_target) < z_tol
order = np.argsort(x_S[mid])
r_s = x_S[mid][order]
sig_mid = sig[mid][order]

lame = analytical_lame_open(r_s)
therm = analytical_thermal_stress(r_s)
ana = [l + t for l, t in zip(lame, therm)]
sig_ana = np.zeros((len(r_s), 3, 3))
for i in range(3):
    sig_ana[:, i, i] = ana[i]


def rel_l2(num, ref):
    return np.sqrt(np.mean((num - ref) ** 2)) / np.sqrt(np.mean(ref**2))


num = [sig_mid[:, i, i] for i in range(3)]
vm_ana, tr_ana = von_mises(sig_ana), tresca(sig_ana)


# --.. ..- .-.. .-.. --- plot --.. ..- .-.. .-.. ---
MPa = 1e-6
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

for i, (lab, col) in enumerate(zip([r"\sigma_{rr}", r"\sigma_{\theta\theta}", r"\sigma_{zz}"],
                                   ["#D55E00", "#009E73", "#0072B2"])):
    ax1.plot(r_s, num[i] * MPa, "o", color=col, ms=4, alpha=0.6, label=rf"Num. ${lab}$")
    ax1.plot(r_s, ana[i] * MPa, "-", color=col, lw=1.5, label=rf"Ana. ${lab}$")
ax1.plot(r_s, np.trace(sig_mid, axis1=1, axis2=2) / 3 * MPa, "k--", lw=1.5, label=r"$p = \mathrm{tr}\,\sigma/3$")
ax1.set_xlabel("Radius (m)")
ax1.set_ylabel("Stress (MPa)")
ax1.set_title("Total stress: Lame + thermal")
ax1.grid(True, ls="--", alpha=0.7)
ax1.legend(frameon=True, fontsize=9)

ax2.plot(r_s, vm_ana * MPa, "-", color="#D55E00", lw=1.5, label=r"Ana. $\sigma_\mathrm{vM}$")
ax2.plot(r_s, von_mises(sig_mid) * MPa, "o", color="#D55E00", ms=4, alpha=0.6, label=r"$\sigma_\mathrm{vM}(\sigma)$")
ax2.plot(r_s, von_mises(sig_mid - np.trace(sig_mid, axis1=1, axis2=2)[:, None, None] / 3 * np.eye(3)) * MPa,
         "x", color="#D55E00", ms=5, label=r"$\sigma_\mathrm{vM}(s)$")
ax2.plot(r_s, tr_ana * MPa, "-", color="#0072B2", lw=1.5, label=r"Ana. $\sigma_\mathrm{T}$")
ax2.plot(r_s, tresca(sig_mid) * MPa, "o", color="#0072B2", ms=4, alpha=0.6, label=r"$\sigma_\mathrm{T}(\sigma)$")
ax2.plot(r_s, tresca(sig_mid - np.trace(sig_mid, axis1=1, axis2=2)[:, None, None] / 3 * np.eye(3)) * MPa,
         "x", color="#0072B2", ms=5, label=r"$\sigma_\mathrm{T}(s)$")
ax2.axhline(0.0, color="k", lw=1.0, ls="--", label=r"$\sigma_\mathrm{vM}(pI) = \sigma_\mathrm{T}(pI) = 0$")
ax2.set_xlabel("Radius (m)")
ax2.set_ylabel("Equivalent stress (MPa)")
ax2.set_title("Equivalent stress of the total and deviatoric parts")
ax2.grid(True, ls="--", alpha=0.7)
ax2.legend(frameon=True, fontsize=9)

fig.suptitle(rf"$P_i$ = {Pi * MPa:.0f} MPa, $T_i$ = {Ti:.0f} K, $T_o$ = {To:.0f} K, open-ended")
fig.tight_layout()
plot_path = os.path.join(CASE_DIR, "output", "invariants.png")
fig.savefig(plot_path, dpi=300)
print(f"[INFO] Plot saved in: {plot_path}")


# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
errors = {
    "L2_error_sigma_rr": error_metric(rel_l2(num[0], ana[0])),
    "L2_error_sigma_tt": error_metric(rel_l2(num[1], ana[1])),
    "L2_error_sigma_zz": error_metric(rel_l2(num[2], ana[2])),
    "L2_error_von_mises": error_metric(rel_l2(von_mises(sig_mid), vm_ana)),
    "L2_error_tresca": error_metric(rel_l2(tresca(sig_mid), tr_ana)),
    "max_von_mises": metric(np.max(von_mises(sig_mid)), np.max(vm_ana)),
    "max_tresca": metric(np.max(tresca(sig_mid)), np.max(tr_ana)),
}
errors.update({k: error_metric(v) for k, v in identities.items()})

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)
