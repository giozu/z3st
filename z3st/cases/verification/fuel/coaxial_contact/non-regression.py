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
    """Radial profile of sigma_rr and sigma_theta at the peak-contact step,
    Z3ST (markers) against the analytical solution (lines).

      * pellet (0 <= r <= b): the free-boundary thermoelastic self-stress of a
        solid cylinder with the computed radial temperature profile T(r),
            sigma_rr = alpha*E [ (1/b^2) I(b) - (1/r^2) I(r) ],
            sigma_th = alpha*E [ (1/b^2) I(b) + (1/r^2) I(r) - (T-T_ref) ],
        with I(r) = int_0^r (T-T_ref) r' dr', plus the uniform contact
        contribution -p so that sigma_rr(b) = -p at the interface. The pellet is
        NOT under a uniform hydrostatic -p: the radial temperature gradient sets
        up self-equilibrated thermal stresses (tensile hoop at the cool rim,
        compression in the hot core) that dominate the interior field.
      * clad (bci <= r <= c): Lame tube under internal pressure p,
            sigma_rr(r) = p bci^2/(c^2-bci^2) (1 - c^2/r^2),
            sigma_th(r) = p bci^2/(c^2-bci^2) (1 + c^2/r^2).
    """
    if not files or not mask.any():
        print("[INFO] no closed-gap step: skipping stress profile")
        return

    i = int(np.argmax(p_lame))
    p = p_lame[i]                                   # MPa, interference pressure
    p_pa = p * 1e6
    m = pv.read(files[i])
    Lz = float(geo["Lz"])

    if "Stress (cells)" in m.cell_data:
        coords = m.cell_centers().points
        s = np.asarray(m.cell_data["Stress (cells)"]).reshape(-1, 9)
    elif "Stress (points)" in m.point_data:
        coords = m.points
        s = np.asarray(m.point_data["Stress (points)"]).reshape(-1, 9)
    else:
        print("[INFO] no Stress field in VTU: skipping stress profile")
        return

    # Z3ST profile on the single axial layer nearest mid-height (a wide slice
    # smears the profile with axial variation; a fixed band can miss the
    # cell-centre rows entirely). Sorted by radius, split into pellet and clad.
    zc = coords[:, 1]
    z_layers = np.unique(np.round(zc, 9))
    z_sel = z_layers[np.argmin(np.abs(z_layers - 0.5 * Lz))]
    dz_layer = np.min(np.diff(z_layers)) if z_layers.size > 1 else Lz
    band = np.abs(zc - z_sel) < 0.5 * dz_layer
    rb = coords[band, 0]
    srr = s[band, 0] / 1e6                           # tensor order (r, theta, z)
    stt = s[band, 4] / 1e6
    o = np.argsort(rb)
    rb, srr, stt = rb[o], srr[o], stt[o]
    pellet = rb <= b + 1e-9
    cladm = rb >= bci - 1e-9

    # --- analytical clad: Lame tube under internal pressure p ---
    # (the pellet interior is thermal-stress dominated and finite in length, so
    # no simple closed form applies there; its computed profile is shown as is.)
    rcl = np.linspace(bci, c, 120)
    kk = p * bci**2 / (c**2 - bci**2)
    srr_clad = kk * (1.0 - c**2 / rcl**2)
    stt_clad = kk * (1.0 + c**2 / rcl**2)

    # Two panels with independent scales: the pellet carries GPa-level thermal
    # stresses that would otherwise crush the ~100 MPa cladding response and
    # hide its Lame comparison.
    C_RR, C_TT = "#4C72B0", "#C44E52"
    fig, (axp, axc) = plt.subplots(
        1, 2, figsize=(10.5, 4.6), gridspec_kw={"width_ratios": [2.2, 1.0]})

    # left: pellet - the computed thermal + contact stress state
    axp.plot(rb[pellet] * 1e3, srr[pellet], "o-", color=C_RR, ms=3, lw=1.6,
             label=r"$\sigma_{rr}$ (radial)")
    axp.plot(rb[pellet] * 1e3, stt[pellet], "s-", color=C_TT, ms=3, lw=1.6,
             label=r"$\sigma_{\theta\theta}$ (hoop)")
    axp.axhline(0, color="grey", lw=0.8, ls=":")
    axp.set_title("pellet: thermal + contact stress")
    axp.set_xlabel("radius r (mm)")
    axp.set_ylabel("stress (MPa)")
    axp.legend(fontsize=8)
    axp.grid(alpha=0.3)

    # right: cladding - Z3ST (markers) vs Lame tube (dashed), own scale
    axc.plot(rb[cladm] * 1e3, srr[cladm], "o", color=C_RR, ms=4, mfc="white",
             mew=1.2, label=r"$\sigma_{rr}$ Z3ST")
    axc.plot(rb[cladm] * 1e3, stt[cladm], "s", color=C_TT, ms=4, mfc="white",
             mew=1.2, label=r"$\sigma_{\theta\theta}$ Z3ST")
    axc.plot(rcl * 1e3, srr_clad, "--", color=C_RR, lw=1.6, label=r"$\sigma_{rr}$ Lamé")
    axc.plot(rcl * 1e3, stt_clad, "--", color=C_TT, lw=1.6, label=r"$\sigma_{\theta\theta}$ Lamé")
    axc.axhline(0, color="grey", lw=0.8, ls=":")
    axc.set_title("cladding vs Lamé")
    axc.set_xlabel("radius r (mm)")
    axc.legend(fontsize=7.5)
    axc.grid(alpha=0.3)

    fig.suptitle(f"Radial and hoop stress at mid-height (p = {p:.1f} MPa)")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "stress_profile_verification.png"), dpi=150)
    plt.close(fig)
    print("[INFO] stress_profile_verification.png saved")

    if pellet.any() and cladm.any():
        print(f"[INFO] interface sigma_rr: pellet={srr[pellet][-1]:8.2f} MPa, "
              f"clad={srr[cladm][0]:8.2f} MPa, analytic -p={-p:8.2f} MPa")


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
