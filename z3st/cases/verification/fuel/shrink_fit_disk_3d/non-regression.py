#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Verification of the penalty contact pressure against the analytical Lame
interference-fit pressure, on a 3D quarter-disk segment (the r-theta quarter
disk extruded along z, with the z0 face as axial mid-plane symmetry and the
z1 face traction-free). The slice is axially free, so — unlike the
plane-strain 2D disk case — the exact solution is the open-ended
(plane-stress, sigma_zz = 0) Lame state. The inner disk is heated uniformly,
so its free thermal expansion is exactly u(b) = alpha_f (T - T_ref) b and the
outer ring does not expand. The interference

    delta = alpha_f (T_pellet - T_ref) * b - g0

then gives the exact Lame shrink-fit pressure for a solid disk in a ring
(plane-stress form):

    p = delta / { b [ (1/E_c)((c^2+b^2)/(c^2-b^2) + nu_c)
                    + (1/E_f)(1 - nu_f) ] }

"""

import os
import glob
import yaml
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv

from z3st.utils.utils_verification import pass_fail_check, regression_check

CASE = os.path.dirname(__file__)
OUT = os.path.join(CASE, "output")
OUT_JSON = os.path.join(CASE, "output", "non-regression.json")

TOLERANCE = 5e-2

files = sorted(glob.glob(os.path.join(OUT, "fields_*.vtu")))

geo = yaml.safe_load(open(os.path.join(CASE, "geometry.yaml")))
b = float(geo["outer_radius_1"])     # pellet outer / interface radius
bci = float(geo["inner_radius_2"])
c = float(geo["outer_radius_2"])

inp = yaml.safe_load(open(os.path.join(CASE, "input.yaml")))
cc = inp["models"]["contact"]
g0 = float(cc["initial_gap"])
fuel = yaml.safe_load(open(os.path.join(CASE, inp["materials"]["cyl_1"])))
Ef, nuf, af, Trf = (float(fuel[k]) for k in ("E", "nu", "alpha", "T_ref"))
clad = yaml.safe_load(open(os.path.join(CASE, inp["materials"]["cyl_2"])))
Ec, nuc = float(clad["E"]), float(clad["nu"])

# Lame interference-fit compliance (plane stress, axially free): delta = p * b * comp
comp = (1.0 / Ec) * ((c**2 + bci**2) / (c**2 - bci**2) + nuc) + (1.0 / Ef) * (1.0 - nuf)

amean = lambda v, r, lo, hi: (lambda m: np.sum(v[m] * r[m]) / np.sum(r[m]))((r >= lo) & (r <= hi))

T_pellet, p_z3st, p_lame = [], [], []
for f in files:
    m = pv.read(f)
    x, y = m.points[:, 0], m.points[:, 1]
    r = np.hypot(x, y)
    T = m.point_data["Temperature"]

    Tp = amean(T, r, 0.0, b)                                  # uniform pellet temperature
    T_pellet.append(Tp)

    p_z3st.append(float(m.cell_data["ContactPressure"][0]) / 1e6)   # MPa

    delta = af * (Tp - Trf) * b - g0                          # exact interference (axially free)
    p_lame.append((delta / (b * comp) / 1e6) if delta > 0 else 0.0)

T_pellet, p_z3st, p_lame = map(np.array, (T_pellet, p_z3st, p_lame))

plt.figure(figsize=(7, 5))
plt.plot(T_pellet, p_lame, "k--", lw=1.5, label="Analytical Lame interference (exact)")
plt.plot(T_pellet, p_z3st, "r-o", lw=2, label="Z3ST penalty contact")
plt.xlabel("pellet temperature (K)")
plt.ylabel("contact pressure (MPa)")
plt.title("Contact pressure verification: Z3ST vs Lame (3D disk segment, uniform ΔT)")
plt.grid(True, ls=":", alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT, "contact_pressure_verification.png"), dpi=150)
print("[INFO] contact_pressure_verification.png saved")

# Compare only where the contact is firmly established: right at closure the
# pressure is a small fraction of peak and the penalty regularization (finite
# k_pen allows ~p/k_pen penetration) dominates the relative error. The
# threshold is relative to the peak so it tracks the material stiffness.
mask = p_lame > 0.25 * p_lame.max()
if not mask.any():
    # p_lame is identically zero, i.e. the interference never becomes
    # positive, so there is no established contact to compare against the
    # Lame fit.
    raise SystemExit(
        "[FAIL] the interference never becomes positive, so contact never "
        "engages; the case premise (the gap closes under the temperature "
        "ramp) does not hold"
    )
rel = np.abs(p_z3st[mask] - p_lame[mask]) / p_lame[mask]
print(f"[INFO] closed-gap steps: {mask.sum()}, mean rel. error vs Lame = {rel.mean() * 100:.1f}%")
print("[INFO] non-regression completed.\n")


def plot_stress_profiles():
    """Radial profile of sigma_rr and sigma_theta at the last step,
    Z3ST vs the analytical Lame interference-fit solution.

      * pellet (0 <= r <= b), solid disk under external pressure p ->
        sigma_rr = sigma_theta = -p   (uniform)
      * clad (bci <= r <= c), ring with internal pressure p, free outer ->
            sigma_rr(r)    = p bci^2/(c^2-bci^2) (1 - c^2/r^2)
            sigma_theta(r) = p bci^2/(c^2-bci^2) (1 + c^2/r^2)

    sigma_rr is continuous (-p) across the gap interface. The Cartesian
    stress tensor is rotated to polar components at each cell center.
    """
    if not files or not mask.any():
        print("[INFO] no closed-gap step: skipping stress profile")
        return

    i = int(np.argmax(p_lame))
    p = p_lame[i]                                   # MPa, exact interference pressure
    m = pv.read(files[i])

    if "Stress (cells)" in m.cell_data:
        coords = m.cell_centers().points
        s = np.asarray(m.cell_data["Stress (cells)"]).reshape(-1, 9)
    elif "Stress (points)" in m.point_data:
        coords = m.points
        s = np.asarray(m.point_data["Stress (points)"]).reshape(-1, 9)
    else:
        print("[INFO] no Stress field in VTU: skipping stress profile")
        return

    xc, yc = coords[:, 0], coords[:, 1]
    rb = np.hypot(xc, yc)
    ct = np.divide(xc, rb, out=np.ones_like(rb), where=rb > 1e-12)
    st = np.divide(yc, rb, out=np.zeros_like(rb), where=rb > 1e-12)
    sxx, syy, sxy = s[:, 0], s[:, 4], s[:, 1]        # Cartesian (x, y, z) tensor
    srr = (sxx * ct**2 + syy * st**2 + 2.0 * sxy * ct * st) / 1e6
    stt = (sxx * st**2 + syy * ct**2 - 2.0 * sxy * ct * st) / 1e6
    o = np.argsort(rb)
    rb, srr, stt = rb[o], srr[o], stt[o]
    pellet = rb <= b + 1e-9
    cladm = rb >= bci - 1e-9

    # analytical curves at the same pressure p
    rp = np.linspace(0.0, b, 50)
    rcl = np.linspace(bci, c, 80)
    kk = p * bci**2 / (c**2 - bci**2)
    srr_clad = kk * (1.0 - c**2 / rcl**2)
    stt_clad = kk * (1.0 + c**2 / rcl**2)

    fig, ax = plt.subplots(figsize=(7.5, 5))
    # analytic (lines)
    ax.plot(rp * 1e3, np.full_like(rp, -p), color="C0", ls="--", lw=1.5,
            label=r"$\sigma_{rr}$ analytic")
    ax.plot(rcl * 1e3, srr_clad, color="C0", ls="--", lw=1.5)
    ax.plot(rp * 1e3, np.full_like(rp, -p), color="C3", ls=":", lw=1.8,
            label=r"$\sigma_{\theta\theta}$ analytic")
    ax.plot(rcl * 1e3, stt_clad, color="C3", ls=":", lw=1.8)
    # Z3ST (markers)
    ax.plot(rb[pellet] * 1e3, srr[pellet], "o", color="C0", ms=4,
            label=r"$\sigma_{rr}$ Z3ST")
    ax.plot(rb[cladm] * 1e3, srr[cladm], "o", color="C0", ms=4)
    ax.plot(rb[pellet] * 1e3, stt[pellet], "s", color="C3", ms=4,
            label=r"$\sigma_{\theta\theta}$ Z3ST")
    ax.plot(rb[cladm] * 1e3, stt[cladm], "s", color="C3", ms=4)

    ax.axvspan(b * 1e3, bci * 1e3, color="0.85", alpha=0.7)
    ax.axhline(0, color="grey", lw=0.8, ls="-")
    ax.text((b + bci) / 2 * 1e3, ax.get_ylim()[1] * 0.9, "gap",
            ha="center", fontsize=8, color="0.4")
    ax.set_xlabel("radius r (mm)")
    ax.set_ylabel("stress (MPa)")
    ax.set_title(f"Radial / hoop stress vs Lame (p = {p:.1f} MPa, 3D disk segment)")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(True, ls=":", alpha=0.5)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "stress_profile_verification.png"), dpi=150)
    plt.close(fig)
    print("[INFO] stress_profile_verification.png saved")

    # interface continuity diagnostic
    if pellet.any() and cladm.any():
        srr_fuel_surf = srr[pellet][-1]
        srr_clad_inner = srr[cladm][0]
        print(f"[INFO] interface sigma_rr: pellet={srr_fuel_surf:8.2f} MPa, "
              f"clad={srr_clad_inner:8.2f} MPa, analytic -p={-p:8.2f} MPa")
        print(f"[INFO] pellet sigma_rr range: [{srr[pellet].min():.2f}, "
              f"{srr[pellet].max():.2f}] MPa (should be ~ -p, uniform, never tensile)")


plot_stress_profiles()

# --. numerical results --..
errors = {
    "contact_pressure": {
        "numerical": p_z3st[mask].max(),
        "reference": p_lame[mask].max(),
        "abs_error": float(np.abs(p_z3st[mask] - p_lame[mask]).max()),   # MPa
        "rel_error": float(rel.max()),
    },
}

pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE)
regression_check(errors, CASE)

print("\n[INFO] non-regression completed.\n")
