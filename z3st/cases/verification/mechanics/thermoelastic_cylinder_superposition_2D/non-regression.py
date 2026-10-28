#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: verification/mechanics/thermoelastic_cylinder_superposition_2D

Open-ended thick cylinder, internal pressure plus a radial temperature
gradient, solved three times on the same mesh:

  total/       Pi = 50 MPa, T = 500/400 K
  thermal/     Pi = 0,      T = 500/400 K    (thermal stresses only)
  mechanical/  Pi = 50 MPa, T = T_ref        (mechanical stresses only)

and twice more with a hydrostatic mechanical load, p = 50 MPa on the inner,
outer and top faces, so that sigma_mech = -p I everywhere:

  total_hydrostatic/       p on every face, T = 500/400 K
  mechanical_hydrostatic/  p on every face, T = T_ref

Checks:

1. each subproblem against its analytic solution (Lame open end; free-end
   thermal stresses);
2. superposition: sigma_tot = sigma_th + sigma_mech, to solver precision
   (linear elasticity, temperature-independent properties);
3. von Mises and Tresca do NOT superpose. Both are seminorms of the stress
   tensor, so the only thing that holds is the triangle inequality
       |f(sigma_th) - f(sigma_mech)| <= f(sigma_tot) <= f(sigma_th) + f(sigma_mech)
   The tensors are added first, the equivalent stress is taken after;
4. the deviatoric/hydrostatic identities of the companion case
   thermoelastic_cylinder_invariants_2D, on each subproblem;
5. with a hydrostatic mechanical load the triangle inequality becomes an
   equality, f(tot) = f(th) + f(mech) = f(th), for any thermal field: the
   hydrostatic part is invisible to both criteria. It is not the only case
   of equality: von Mises also has it for proportional deviators, Tresca
   wherever both states share principal directions and their ordering
   (tracked below on the Pi case).
"""

import os

import matplotlib.pyplot as plt
import numpy as np

from z3st.utils.non_regression import error_metric, finish, load_yaml, metric, tracked
from z3st.utils.utils_extract_vtu import extract_field

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")
os.makedirs(os.path.join(CASE_DIR, "output"), exist_ok=True)
SUBS = ("total", "thermal", "mechanical", "total_hydrostatic", "mechanical_hydrostatic")
PAIRS = {"total": "mechanical", "total_hydrostatic": "mechanical_hydrostatic"}  # total -> its mechanical part
EQUALITY = 1.0e-8  # -          tolerance for f(tot) = f(th) + f(mech), hydrostatic load

geom = load_yaml(CASE_DIR, "geometry.yaml")
mat = load_yaml(CASE_DIR, os.path.join("total", load_yaml(CASE_DIR, "total/input.yaml")["materials"]["steel"]))
Ri, Ro, Lz = float(geom["Ri"]), float(geom["Ro"]), float(geom["Lz"])  # m   inner/outer radius, height
E, nu, alpha = float(mat["E"]), float(mat["nu"]), float(mat["alpha"])  # Pa, -, 1/K

Pi, Po = 50.0e6, 0.0  # Pa         internal and external pressure
Ti, To = 500.0, 400.0  # K          inner and outer surface temperature
z_target, z_tol = Lz / 2, Lz / 10  # m          z-plane for data extraction

TOLERANCE = 3.0e-2  # -          tolerance for non-regression
ROUNDOFF = 1.0e-12  # -          tolerance for the algebraic identities
SUPERPOSITION = 1.0e-8  # -       tolerance for sigma_tot = sigma_th + sigma_mech


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


def diag(components):
    """(N, 3, 3) tensor from the three diagonal profiles."""
    out = np.zeros((len(components[0]), 3, 3))
    for i in range(3):
        out[:, i, i] = components[i]
    return out


def hydrostatic(sig):
    return np.trace(sig, axis1=-2, axis2=-1)[..., None, None] / 3.0 * np.eye(3)


def von_mises(sig):
    s = sig - hydrostatic(sig)
    return np.sqrt(1.5 * np.einsum("...ij,...ij", s, s))


def tresca(sig):
    """sigma_1 - sigma_3 from the principal stresses."""
    ev = np.linalg.eigvalsh(sig)
    return ev[..., -1] - ev[..., 0]


def rel_l2(num, ref):
    return np.sqrt(np.mean((num - ref) ** 2)) / np.sqrt(np.mean(ref**2))


# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
sig, x_S, z_S = {}, None, None
for sub in SUBS:
    x, z, _, S = extract_field(os.path.join(CASE_DIR, sub, "output", "fields.vtu"), field_name="Stress (cells)")
    if x_S is not None:
        assert np.array_equal(x, x_S) and np.array_equal(z, z_S), f"{sub}: cell order differs from total"
    x_S, z_S = x, z
    sig[sub] = S.reshape(-1, 3, 3)  # axisymmetric layout (r, theta, z)

vm = {k: von_mises(v) for k, v in sig.items()}
tr = {k: tresca(v) for k, v in sig.items()}

# 2. superposition of the tensors, whole domain, for both mechanical loads
superposition = {}
for tot_k, mech_k in PAIRS.items():
    superposition[tot_k] = np.max(np.abs(sig[tot_k] - sig["thermal"] - sig[mech_k])) / np.max(np.abs(sig[tot_k]))
    print(f"[INFO] {tot_k}: max |sigma_tot - sigma_th - sigma_mech| / max|sigma_tot| = {superposition[tot_k]:.3e}")
    assert superposition[tot_k] < SUPERPOSITION, f"{tot_k}: superposition broken"

# The hydrostatic load really is hydrostatic: sigma_mech = -p I
hydro_dev = np.max(np.abs(sig["mechanical_hydrostatic"] + Pi * np.eye(3))) / Pi
print(f"[INFO] max |sigma_mech_hyd + p I| / p = {hydro_dev:.3e}")
assert hydro_dev < SUPERPOSITION, "the hydrostatic mechanical load is not -p I"

# 3. the equivalent stresses only obey the triangle inequality ...
# 5. ... which becomes an equality when the mechanical load is hydrostatic
equality = {}
for name, f in (("von Mises", vm), ("Tresca", tr)):
    th = f["thermal"]
    for tot_k, mech_k in PAIRS.items():
        tot, me = f[tot_k], f[mech_k]
        slack = SUPERPOSITION * np.max(tot)  # superposition is exact only to solver precision
        assert np.all(tot <= th + me + slack), f"{name}, {tot_k}: f(tot) > f(th) + f(mech)"
        assert np.all(tot >= np.abs(th - me) - slack), f"{name}, {tot_k}: f(tot) < |f(th) - f(mech)|"
        gap = (th + me - tot) / np.max(tot)
        print(f"[INFO] {name:9s}, {tot_k:17s}: f(th) + f(mech) - f(tot) up to {np.max(gap) * 100:.3g}% "
              f"of max f(tot)")
    equality[name] = np.max(np.abs(f["total_hydrostatic"] - th - f["mechanical_hydrostatic"])) / np.max(th)
    assert equality[name] < EQUALITY, f"{name}: equality broken under a hydrostatic load"

# Equality without any hydrostatic load: Tresca on the Pi case, where the
# thermal and mechanical states share principal directions and ordering
coaxial = np.abs(tr["thermal"] + tr["mechanical"] - tr["total"]) < EQUALITY * np.max(tr["total"])
print(f"[INFO] Tresca, Pi case: f(tot) = f(th) + f(mech) in {np.mean(coaxial) * 100:.1f}% of the cells")

# 4. deviatoric / hydrostatic identities on every subproblem
identities = {}
for sub in SUBS:
    s, h = sig[sub] - hydrostatic(sig[sub]), hydrostatic(sig[sub])
    ref = np.max(vm["total"])  # not vm[sub]: the hydrostatic load has vM ~ 0
    identities[f"{sub}_vM_dev_minus_vM_tot"] = np.max(np.abs(von_mises(s) - vm[sub])) / ref
    identities[f"{sub}_vM_hyd"] = np.max(von_mises(h)) / ref
    identities[f"{sub}_Tresca_dev_minus_Tresca_tot"] = np.max(np.abs(tresca(s) - tr[sub])) / ref
    identities[f"{sub}_Tresca_hyd"] = np.max(tresca(h)) / ref
for k, v in identities.items():
    print(f"[INFO] {k:40s} = {v:.3e}")
    assert v < ROUNDOFF, f"{k} = {v:.3e} exceeds round-off ({ROUNDOFF:.0e})"

# 1. mid-plane profiles against the analytic solutions
mid = np.abs(z_S - z_target) < z_tol
order = np.argsort(x_S[mid])
r_s = x_S[mid][order]
mid_sig = {k: v[mid][order] for k, v in sig.items()}

therm, lame = analytical_thermal_stress(r_s), analytical_lame_open(r_s)
ana = {
    "thermal": diag(therm),
    "mechanical": diag(lame),
    "total": diag([t + m for t, m in zip(therm, lame)]),
    "mechanical_hydrostatic": diag([np.full_like(r_s, -Pi)] * 3),
    "total_hydrostatic": diag([t - Pi for t in therm]),
}


# --.. ..- .-.. .-.. --- plot --.. ..- .-.. .-.. ---
MPa = 1e-6
colours = {"total": "k", "thermal": "#D55E00", "mechanical": "#0072B2"}
fig, axes = plt.subplots(2, 2, figsize=(13, 10), sharex=True)

rows = (("total", "mechanical", rf"$P_i$ = {Pi * MPa:.0f} MPa on the inner face"),
        ("total_hydrostatic", "mechanical_hydrostatic", rf"$p$ = {Pi * MPa:.0f} MPa on every face ($\sigma = -pI$)"))
for row, (tot_k, mech_k, load) in zip(axes, rows):
    for ax, name, f in zip(row, ("von Mises", "Tresca"), (von_mises, tresca)):
        for sub, key in (("total", tot_k), ("thermal", "thermal"), ("mechanical", mech_k)):
            ax.plot(r_s, f(ana[key]) * MPa, "-", color=colours[sub], lw=1.5)
            ax.plot(r_s, f(mid_sig[key]) * MPa, "o", color=colours[sub], ms=4, alpha=0.6, label=sub)
        th, me = f(mid_sig["thermal"]), f(mid_sig[mech_k])
        ax.plot(r_s, (th + me) * MPa, "--", color="#009E73", lw=2, label="thermal + mechanical")
        ax.plot(r_s, np.abs(th - me) * MPa, ":", color="#009E73", lw=2, label="|thermal - mechanical|")
        ax.set_title(f"{name}, {load}", fontsize=11)
        ax.grid(True, ls="--", alpha=0.7)
        ax.legend(frameon=True, fontsize=9)
    row[0].set_ylabel("Equivalent stress (MPa)")
for ax in axes[1]:
    ax.set_xlabel("Radius (m)")

fig.suptitle(rf"$T_i$ = {Ti:.0f} K, $T_o$ = {To:.0f} K; lines analytic, markers FE. "
             "Top: the sum overestimates. Bottom: hydrostatic load, the sum is exact")
fig.tight_layout()
plot_path = os.path.join(CASE_DIR, "output", "superposition.png")
fig.savefig(plot_path, dpi=300)
print(f"[INFO] Plot saved in: {plot_path}")


# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
errors = {}
for sub in ("total", "thermal", "mechanical", "total_hydrostatic"):
    for i, c in enumerate(("rr", "tt", "zz")):
        if np.any(ana[sub][:, i, i]):  # sigma_zz of the mechanical subproblem is 0 analytically
            errors[f"{sub}_L2_error_sigma_{c}"] = error_metric(rel_l2(mid_sig[sub][:, i, i], ana[sub][:, i, i]))
    errors[f"{sub}_L2_error_von_mises"] = error_metric(rel_l2(von_mises(mid_sig[sub]), von_mises(ana[sub])))
    errors[f"{sub}_L2_error_tresca"] = error_metric(rel_l2(tresca(mid_sig[sub]), tresca(ana[sub])))
    errors[f"{sub}_max_von_mises"] = metric(np.max(von_mises(mid_sig[sub])), np.max(von_mises(ana[sub])))
errors["mechanical_max_abs_sigma_zz_over_Pi"] = error_metric(np.max(np.abs(mid_sig["mechanical"][:, 2, 2])) / Pi)
errors.update({f"superposition_{k}": error_metric(v) for k, v in superposition.items()})
errors["mechanical_hydrostatic_max_dev_from_minus_pI"] = error_metric(hydro_dev)
errors.update({f"equality_hydrostatic_{k.replace(' ', '_')}": error_metric(v) for k, v in equality.items()})
errors["Tresca_coaxial_fraction_Pi_case"] = tracked(np.mean(coaxial))
errors["vM_sum_overestimate_max"] = tracked(np.max((vm["thermal"] + vm["mechanical"] - vm["total"]) / np.max(vm["total"])))
errors.update({k: error_metric(v) for k, v in identities.items()})

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)
