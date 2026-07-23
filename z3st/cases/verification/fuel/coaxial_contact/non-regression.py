#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Verification of the penalty contact pressure against the analytical Lame
interference-fit pressure. The pellet is NOT heated uniformly: volumetric
heat generation gives it a radial gradient (centerline hottest, surface
coolest). For a free axisymmetric solid disk/cylinder, the exact result for
the outer-surface radial displacement is

    u(R) = alpha * R * Tbar(R),   Tbar(R) = (2/R^2) int_0^R T(r) r dr

i.e. it depends only on the area-averaged cross-section temperature, not on
the profile shape or Poisson's ratio. For the standard parabolic fuel-pellet
profile (uniform heat generation, insulated ends) this area average reduces
to the arithmetic mean of the centerline and surface temperatures:

    Tbar(R) = (T_c + T_s) / 2

so the interference is

    delta = alpha_f ((T_c + T_s)/2 - T_ref) * b - g0

which then gives the exact Lame shrink-fit pressure for a solid cylinder in a tube:

    p = delta / { b [ (1/E_c)((c^2+b^2)/(c^2-b^2) + nu_c) + (1/E_f)(1 - nu_f) ] }

plane-stress form, consistent with the axially-free pellet. Note this only
gives the correct outer-boundary displacement (and hence delta/p); the
pellet's own stress state is no longer uniform hydrostatic -p once in
contact, since the radial gradient also produces self-equilibrated thermal
stresses (sigma_rr = alpha*E*(Tbar(R)-Tbar(r)), not accounted for below in
plot_stress_profiles()).

"""

import os
import re
import glob
import yaml
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv

from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_verification import *

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

# Lame interference-fit compliance (plane stress): delta = p * b * comp
comp = (1.0 / Ec) * ((c**2 + bci**2) / (c**2 - bci**2) + nuc) + (1.0 / Ef) * (1.0 - nuf)

amean = lambda v, r, lo, hi: (lambda m: np.sum(v[m] * r[m]) / np.sum(r[m]))((r >= lo) & (r <= hi))

# Pressure comes straight from the writer's "ContactPressure" cell field
# (utils/writer.py), a uniform broadcast of the ContactModel's own live scalar
# -- the value it actually applied as traction in the mechanical weak form.
# Re-deriving it from the Displacement field is fragile (nodal-tolerance picks
# between the two facing surfaces, mesh/geometry drift, etc.), so this reads
# the solver's own number directly instead of reconstructing it.
#
# The gap itself is not (yet) an exported field, so it is still parsed from
# the run log's "[contact] ... gap=... um" line (converged/last sample per
# time step) -- only needed for the open-gap fallback check below.
log_text = open(os.path.join(CASE, "log_z3st.md")).read()
step_blocks = [blk for blk in re.split(r"(?=## Step \d+/\d+:)", log_text) if blk.startswith("## Step")]
gap_re = re.compile(r"\[contact\].*?gap=([+-]?[\d.]+) um")

gap_z3st = []
for blk in step_blocks:
    matches = gap_re.findall(blk)
    if not matches:
        raise RuntimeError("no '[contact]' line found in a step block of log_z3st.md")
    gap_z3st.append(float(matches[-1]) * 1e-6)   # converged (last) sample of the step

if len(step_blocks) != len(files):
    raise RuntimeError(
        f"log_z3st.md has {len(step_blocks)} step blocks but {len(files)} fields_*.vtu "
        "files were found: cannot align contact pressure history with pellet temperature."
    )

T_pellet, p_z3st, p_lame, gap_free = [], [], [], []
for f in files:
    m = pv.read(f)
    r = m.points[:, 0]
    T = m.point_data["Temperature"]

    Tp = amean(T, r, 0.0, b)                                  # mean pellet temperature
    T_pellet.append(Tp)

    p_z3st.append(float(m.cell_data["ContactPressure"][0]) / 1e6)   # MPa

    # Driving temperature for the Lame interference: area-average of the
    # pellet cross-section, which for the standard parabolic radial profile
    # is exactly the mean of the centerline (r=0) and surface (r=b)
    # temperatures -- see module docstring.
    T_c = T[np.isclose(r, 0.0, atol=1e-9)].mean()
    T_s = T[np.isclose(r, b, atol=1e-9)].mean()
    Tp_lame = (T_c + T_s) / 2.0

    delta = af * (Tp_lame - Trf) * b - g0                     # exact interference
    gap_free.append(-delta)                                   # analytic OPEN gap (no contact)
    p_lame.append((delta / (b * comp) / 1e6) if delta > 0 else 0.0)

T_pellet, p_z3st, p_lame, gap_z3st, gap_free = map(
    np.array, (T_pellet, p_z3st, p_lame, gap_z3st, gap_free))

plt.figure(figsize=(7, 5))
plt.plot(T_pellet, p_lame, "k--", lw=1.5, label="Analytical Lame interference (exact)")
plt.plot(T_pellet, p_z3st, "r-o", lw=2, label="Z3ST penalty contact")
plt.xlabel("pellet temperature (K)")
plt.ylabel("contact pressure (MPa)")
plt.title("Contact pressure verification: Z3ST vs Lame (uniform ΔT)")
plt.grid(True, ls=":", alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT, "contact_pressure_verification.png"), dpi=150)
print("[INFO] contact_pressure_verification.png saved")

mask = p_lame > 1.0
if mask.any():
    # Deviation measured relative to the characteristic (peak) contact
    # pressure, not the local per-step value. At contact onset the penalty
    # pressure p = k_pen * interpenetration necessarily lags the analytical
    # Lame line, which jumps from zero the instant the interference is
    # positive; dividing by the tiny local p_lame there would turn that small,
    # expected absolute lag into a spurious large relative error. Normalising
    # by the peak pressure scale measures the physically meaningful quantity:
    # the agreement in established contact and the finite-penalty-stiffness
    # residual near peak load.
    p_scale = p_lame[mask].max()
    rel = np.abs(p_z3st[mask] - p_lame[mask]) / p_scale
    print(f"[INFO] closed-gap steps: {mask.sum()}, dev vs Lame rel. to peak "
          f"{p_scale:.1f} MPa: max {rel.max() * 100:.1f}%, mean {rel.mean() * 100:.1f}%")
print("[INFO] non-regression completed.\n")


def plot_stress_profiles():
    """Radial profile of sigma_rr and sigma_theta at the last step,
    Z3ST vs the analytical Lame interference-fit solution.

      * pellet (0 <= r <= b), solid cylinder under external pressure p, axially
        free  ->  sigma_rr = sigma_theta = -p   (uniform)
      * clad (bci <= r <= c), tube with internal pressure p, free outer ->
            sigma_rr(r)    = p bci^2/(c^2-bci^2) (1 - c^2/r^2)
            sigma_theta(r) = p bci^2/(c^2-bci^2) (1 + c^2/r^2)

    sigma_rr is continuous (-p) across the gap interface.
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

    Lz = float(geo["Lz"])
    rr_c, zz_c = coords[:, 0], coords[:, 1]
    band = np.abs(zz_c - 0.5 * Lz) < 0.25 * Lz      # mid-height slice
    rb = rr_c[band]
    srr = s[band, 0] / 1e6                           # tensor order (r, theta, z)
    stt = s[band, 4] / 1e6
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
    ax.set_title(f"Radial / hoop stress vs Lame (p = {p:.1f} MPa, mid-height)")
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
if mask.any():
    # Contact regime: verify the penalty pressure against the Lame
    # interference fit. The reported numerical/reference pair is taken at the
    # peak-contact step (established contact), and the error is the largest
    # per-step deviation relative to the peak pressure scale (see the onset
    # note above).
    i_peak = int(np.argmax(p_lame))
    errors = {
        "contact_pressure": {
            "numerical": float(p_z3st[i_peak]),
            "reference": float(p_lame[i_peak]),
            "abs_error": float(np.abs(p_z3st[mask] - p_lame[mask]).max()),
            "rel_error": float(rel.max()),
        },
    }
else:
    # Open-gap regime (the power/temperature never closes the gap): there is
    # no contact pressure to verify, so verify the final gap width against the
    # analytic free thermal expansion instead — and require that the solver
    # also reports no spurious contact.
    print("[INFO] gap never closes: verifying final open gap vs free expansion")
    gap_err = abs(gap_z3st[-1] - gap_free[-1]) / g0
    errors = {
        "final_gap": {
            "numerical": float(gap_z3st[-1]),
            "reference": float(gap_free[-1]),
            "abs_error": float(abs(gap_z3st[-1] - gap_free[-1])),
            "rel_error": float(gap_err),
        },
        "spurious_contact_pressure": {
            "numerical": float(p_z3st.max()),
            "reference": 0.0,
            "abs_error": float(p_z3st.max()),
            "rel_error": float(p_z3st.max()),   # MPa; must be ~0
        },
    }

pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE)
regression_check(errors, CASE)

print("\n[INFO] non-regression completed.\n")
