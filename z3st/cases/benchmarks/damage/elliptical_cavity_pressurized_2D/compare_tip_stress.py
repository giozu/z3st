#!/usr/bin/env python3
"""
Compare the local stress state ahead of the cavity tip, at the moment of
crack initiation (max(Damage) >= DAMAGE_THRESHOLD), between:

  - elliptical_cavity_2D              (remote displacement loading, ymax)
  - elliptical_cavity_pressurized_2D  (internal pressure loading, cavity)

Both cases now share the same lenticular cavity geometry (ax=11.2um,
ay=5um) and the same calibrated material (uo2.yaml, sigma_c=150 MPa), so
this isolates the effect of the loading mode alone: does the same local
stress trigger fracture in both, independent of how the load is applied?

The two raw "headline" numbers (188 MPa remote reaction stress / 177 MPa
internal cavity pressure) are not directly comparable -- this script
extracts the actual local stress ahead of the cavity tip in both cases so
they can be compared on equal footing against sigma_c.

Where does cracking actually start? It turns out NOT at the same point in
both cases: the original (remote tension) initiates exactly at the
major-axis tip (ax, 0) as classical Inglis theory predicts, but the copy
(internal pressure) initiates at a point offset ~28 deg around the
boundary -- pure internal pressurization of an elongated cavity is known
to shift the tensile concentration away from the major-axis tip, unlike
remote tension. So each case's stress profile below is anchored on ITS
OWN actual max(Damage) location (not an assumed point), and probed
radially outward from the cavity center through that point.

That anchor point can also sit on (or very near) a geometric corner --
the two circular arcs of MakeLenticularbubble (mesh.geo) do NOT meet
tangentially at (ax, 0) (their tangent directions differ by ~12.5 deg) --
so raw nodal stress right at such a point can be singular (Williams'
corner solution) and mesh-dependent. Hence the stress-vs-distance sweep
(multiples of lc, the model's own regularization length) rather than a
single-point reading -- the standard "Point/Line Method" for notch stress
fields.
"""

import os
import re

import numpy as np
import yaml
import matplotlib.pyplot as plt

from z3st.utils.utils_extract_vtu import extract_field

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
ORIGINAL_DIR = os.path.normpath(os.path.join(THIS_DIR, "..", "elliptical_cavity_2D"))
COPY_DIR = THIS_DIR

DAMAGE_THRESHOLD = 0.5
SIGMA_C_MPA = 150.0
STANDOFF_MULTIPLES_OF_LC = [0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0]


def step_of(fname):
    m = re.search(r'_(\d+)\.vtu', fname)
    return int(m.group(1)) if m else -1


def find_critical_vtu(case_dir):
    """First simulation_*.vtu where max(Damage) over the whole domain
    crosses DAMAGE_THRESHOLD; falls back to the last available step if it
    never crosses (with a warning)."""
    out_dir = os.path.join(case_dir, "output")
    vtu_files = sorted(
        (f for f in os.listdir(out_dir) if f.startswith("simulation_") and f.endswith(".vtu")),
        key=step_of,
    )
    if not vtu_files:
        raise FileNotFoundError(f"No simulation_*.vtu in {out_dir}")

    for fname in vtu_files:
        path = os.path.join(out_dir, fname)
        _, _, _, D = extract_field(path, field_name="Damage")
        if np.max(D) >= DAMAGE_THRESHOLD:
            return path, step_of(fname)

    print(f"[WARNING] {case_dir}: max(Damage) never reached {DAMAGE_THRESHOLD}; "
          f"using the last available step instead.")
    path = os.path.join(out_dir, vtu_files[-1])
    return path, step_of(vtu_files[-1])


def von_mises(S):
    """S: (N, 9) row-major flattened 3x3 stress tensor per point."""
    sxx, sxy, sxz = S[:, 0], S[:, 1], S[:, 2]
    syy, syz = S[:, 4], S[:, 5]
    szz = S[:, 8]
    return np.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2
                           + 6.0 * (sxy ** 2 + syz ** 2 + sxz ** 2)))


def read_lc(case_dir):
    with open(os.path.join(case_dir, "input.yaml")) as f:
        cfg = yaml.safe_load(f)
    return float(cfg["damage"]["lc"])


def applied_load_at_step(case_dir, step):
    """Read boundary_conditions.yaml and return (label, value) of the
    ramped load at the given step index, whichever BC type this case uses."""
    with open(os.path.join(case_dir, "boundary_conditions.yaml")) as f:
        bc = yaml.safe_load(f)
    for entry in bc["mechanical"]["uo2"]:
        if entry.get("type") == "Neumann":
            return "cavity pressure (MPa)", float(entry["traction"][step]) / 1e6
        if entry.get("type") == "Dirichlet_y":
            return "ymax displacement (nm)", float(entry["value"][step]) * 1e9
    return "load", float("nan")


def probe(x, y, S, xy):
    """Stress components at the mesh node nearest to point xy."""
    dist = np.hypot(x - xy[0], y - xy[1])
    i = np.argmin(dist)
    return {
        "dist_to_target_nm": float(dist[i]) * 1e9,
        "sigma_xx": float(S[i, 0]),
        "sigma_yy": float(S[i, 4]),
        "von_mises": float(von_mises(S[i:i + 1])[0]),
    }


def analyze_case(case_dir, lc, ax, ay):
    """Find the actual max(Damage) location at the critical step, then
    probe the stress field radially outward from the cavity center through
    that point, at increasing standoff distances (multiples of lc)."""
    vtu_path, step = find_critical_vtu(case_dir)

    x, y, _, D = extract_field(vtu_path, field_name="Damage")
    D = D.reshape(-1)
    i_max = np.argmax(D)
    anchor = np.array([x[i_max], y[i_max]])
    d_max = float(D[i_max])
    theta_deg = float(np.degrees(np.arctan2(anchor[1], anchor[0])))

    direction = anchor / np.linalg.norm(anchor)  # radially outward from (0,0)

    _, _, _, S = extract_field(vtu_path, field_name="Stress_uo2 (points)")
    distances = np.array(STANDOFF_MULTIPLES_OF_LC) * lc
    raw_anchor = probe(x, y, S, tuple(anchor))
    rows = [probe(x, y, S, tuple(anchor + d * direction)) for d in distances]

    load_label, load_value = applied_load_at_step(case_dir, step)

    return {
        "vtu": os.path.basename(vtu_path),
        "step": step,
        "anchor": anchor,
        "theta_deg": theta_deg,
        "d_max": d_max,
        "load_label": load_label,
        "load_value": load_value,
        "raw_anchor": raw_anchor,
        "distances": distances,
        "rows": rows,
        "ax": ax,
        "ay": ay,
    }


def print_profile(label, res):
    print("=" * 78)
    print(f" {label}")
    print("=" * 78)
    print(f"Critical step   : {res['step']} ({res['vtu']})")
    print(f"Applied load    : {res['load_value']:.3f} {res['load_label']}")
    a = res["anchor"]
    print(f"Crack-initiation point (max Damage={res['d_max']:.4f}): "
          f"({a[0]*1e6:.3f}, {a[1]*1e6:.3f}) um, angle {res['theta_deg']:.1f} deg from x-axis")
    print(f"  for reference: major-axis tip (ax,0)=({res['ax']*1e6:.3f}, 0.000) um  |  "
          f"minor-axis tip (0,ay)=(0.000, {res['ay']*1e6:.3f}) um")
    rc = res["raw_anchor"]
    print(f"Raw value AT that point ({rc['dist_to_target_nm']:.2f} nm off) -- "
          f"possibly singular, reference only:")
    print(f"  sigma_xx={rc['sigma_xx']/1e6:.2f} MPa  sigma_yy={rc['sigma_yy']/1e6:.2f} MPa  "
          f"von Mises={rc['von_mises']/1e6:.2f} MPa\n")
    print(f"{'dist (um)':>10}{'dist/lc':>10}{'node off (nm)':>14}"
          f"{'sigma_xx (MPa)':>16}{'sigma_yy (MPa)':>16}{'von Mises (MPa)':>18}")
    for mult, d, row in zip(STANDOFF_MULTIPLES_OF_LC, res["distances"], res["rows"]):
        print(f"{d*1e6:>10.3f}{mult:>10.2f}{row['dist_to_target_nm']:>14.1f}"
              f"{row['sigma_xx']/1e6:>16.2f}{row['sigma_yy']/1e6:>16.2f}{row['von_mises']/1e6:>18.2f}")
    print()


def plot_profiles(distances, res_orig, res_copy, out_path):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    d_um = distances * 1e6

    ax = axes[0]
    ax.plot(d_um, [r["sigma_yy"] / 1e6 for r in res_orig["rows"]], 'o-',
            color='tab:blue', label="Original (displacement)")
    ax.plot(d_um, [r["sigma_yy"] / 1e6 for r in res_copy["rows"]], 'o-',
            color='tab:red', label="Copy (pressure)")
    ax.axhline(SIGMA_C_MPA, color='gray', ls='--', lw=1, label=f"sigma_c = {SIGMA_C_MPA:.0f} MPa")
    ax.axhline(0, color='k', lw=0.7)
    ax.set_xscale('log')
    ax.set_xlabel("Distance from crack-initiation point, radially outward (um)")
    ax.set_ylabel("sigma_yy (MPa)")
    ax.set_title("sigma_yy vs. distance from each case's own crack-initiation point")
    ax.grid(True, ls=':', alpha=0.5, which='both')
    ax.legend(fontsize=8)

    ax = axes[1]
    ax.plot(d_um, [r["von_mises"] / 1e6 for r in res_orig["rows"]], 'o-',
            color='tab:blue', label="Original (displacement)")
    ax.plot(d_um, [r["von_mises"] / 1e6 for r in res_copy["rows"]], 'o-',
            color='tab:red', label="Copy (pressure)")
    ax.axhline(SIGMA_C_MPA, color='gray', ls='--', lw=1, label=f"sigma_c = {SIGMA_C_MPA:.0f} MPa")
    ax.set_xscale('log')
    ax.set_xlabel("Distance from crack-initiation point, radially outward (um)")
    ax.set_ylabel("von Mises (MPa)")
    ax.set_title("von Mises vs. distance from each case's own crack-initiation point")
    ax.grid(True, ls=':', alpha=0.5, which='both')
    ax.legend(fontsize=8)

    fig.suptitle("Stress profile radially outward from each case's own crack-initiation point")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"[INFO] Plot saved: {out_path}")


def main():
    with open(os.path.join(COPY_DIR, "geometry.yaml")) as f:
        geom = yaml.safe_load(f)
    ax = float(geom["cavity"]["ax"])
    ay = float(geom["cavity"]["ay"])

    lc_orig = read_lc(ORIGINAL_DIR)
    lc_copy = read_lc(COPY_DIR)
    if abs(lc_orig - lc_copy) > 1e-15:
        print(f"[WARNING] lc differs between cases ({lc_orig:.3e} vs {lc_copy:.3e}); "
              f"using the copy's lc for both profiles.")
    lc = lc_copy

    print(f"[INFO] Each case is probed at its OWN max(Damage) location, radially "
          f"outward from the cavity center, at {STANDOFF_MULTIPLES_OF_LC} x lc "
          f"(lc = {lc * 1e6:.3f} um).\n")

    res_orig = analyze_case(ORIGINAL_DIR, lc, ax, ay)
    print_profile("original (displacement)", res_orig)

    res_copy = analyze_case(COPY_DIR, lc, ax, ay)
    print_profile("copy (pressure)", res_copy)

    plot_profiles(res_orig["distances"], res_orig, res_copy,
                  os.path.join(COPY_DIR, "ligament_stress_profile.png"))


if __name__ == "__main__":
    main()
