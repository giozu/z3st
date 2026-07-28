#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: verification/fuel/creep_shrink_fit_2D

PCMI coupling chain (thermal + mechanical + Gas gap
conductance + penalty contact + burnup state bus + swelling eigenstrain bus).
All metrics are read from ``output/history.csv`` — the
per-step trajectory streamed by the case-local ``diagnostics.py`` — so the
check is independent of the output format (the run writes a single XDMF).

One metric has a closed form and is checked analytically:

  * ``burnup_avg_final`` — the nodal-mean burnup over the fuel equals the flat
    closed form  bu = Σ_k lhr_k·Δt_k / (area·ρ·HM·8.64e10), because the
    radial form factor is area-normalised to mean 1 (the source bus preserves
    the rating) and ``update_state`` accumulates with the right-endpoint rule
    over the generated power history.

The PCMI end-state scalars (gap, contact pressure, temperatures) have no
closed form — they are recorded with ``rel_error = 0`` so the analytic
pass/fail gate ignores them, and are protected purely by the gold comparison
(``regression_check`` vs ``output/non-regression_gold.json``).
"""

import os
import csv
import yaml
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from z3st.utils.utils_extract_xdmf import extract_field_xdmf, list_fields_xdmf
from z3st.utils.utils_verification import pass_fail_check, regression_check
from z3st.utils.utils_load import generate_power_history

CASE_DIR = os.path.dirname(__file__)
OUT = os.path.join(CASE_DIR, "output")
OUT_JSON = os.path.join(OUT, "non-regression.json")
HISTORY = os.path.join(OUT, "history.csv")
XDMF_FILE = os.path.join(OUT, "fields.xdmf")

TOLERANCE = 1e-2

# --. trajectory from the case-local diagnostics CSV --..
with open(HISTORY) as f:
    rows = list(csv.DictReader(f))
if not rows:
    raise RuntimeError(f"{HISTORY} is empty — did the run complete?")
last = rows[-1]

bu_avg = float(last["burnup_avg_MWdkgU"])
bu_max = float(last["burnup_max_MWdkgU"])
gap_um = float(last["gap_um"])
p_mpa = float(last["contact_pressure_MPa"])
t_max_final = float(last["T_max_K"])
t_max_peak = max(float(r["T_max_K"]) for r in rows)

# --. closed-form mean burnup from the input power history --..
with open(os.path.join(CASE_DIR, "input.yaml")) as f:
    inp = yaml.safe_load(f)
with open(os.path.join(CASE_DIR, "geometry.yaml")) as f:
    geom = yaml.safe_load(f)
fuel_card = inp["materials"]["fuel"]
with open(os.path.join(CASE_DIR, fuel_card)) as f:
    fuel = yaml.safe_load(f)

Ro = float(geom["outer_radius_1"])
area = np.pi * Ro**2
rho = float(fuel["rho"])
hm = float(fuel.get("heavy_metal_fraction", 0.8815))
SECONDS_PER_MWD = 8.64e10  # 86400 s/day * 1e6 W/MW

raw_n_steps = inp["n_steps"]
n_increments = raw_n_steps if isinstance(raw_n_steps, (list, tuple)) else int(raw_n_steps) - 1
times, lhrs, _ = generate_power_history(
    inp["time"], inp["lhr"], n_steps=n_increments, filename=None
)
energy_per_m = float(np.sum(np.asarray(lhrs)[1:] * np.diff(np.asarray(times))))
BU_REF = energy_per_m / (area * rho * hm * SECONDS_PER_MWD)

print(f"[INFO] final mean burnup : numerical = {bu_avg:.4f}, "
      f"closed form = {BU_REF:.4f} MWd/kgU")
print(f"[INFO] final peak burnup : {bu_max:.4f} MWd/kgU (rim)")
print(f"[INFO] final gap         : {gap_um:.4f} um "
      f"({'closed — PCMI active' if gap_um <= 0 else 'open'})")
print(f"[INFO] contact pressure  : {p_mpa:.4f} MPa")
print(f"[INFO] T_max final / peak: {t_max_final:.2f} / {t_max_peak:.2f} K")

if gap_um > 0:
    print("[WARNING] gap did not close by end of run — PCMI path not exercised.")

def _regression_only(value):
    """Metric with no closed form: gold-comparison protection only."""
    return {"numerical": float(value), "reference": None,
            "abs_error": 0.0, "rel_error": 0.0}

errors = {
    "burnup_avg_final": {
        "numerical": bu_avg,
        "reference": BU_REF,
        "abs_error": float(abs(bu_avg - BU_REF)),
        "rel_error": float(abs(bu_avg - BU_REF) / BU_REF),
    },
    "burnup_max_final": _regression_only(bu_max),
    "gap_final_um": _regression_only(gap_um),
    "contact_pressure_final_MPa": _regression_only(p_mpa),
    "T_max_final_K": _regression_only(t_max_final),
    "T_max_peak_K": _regression_only(t_max_peak),
}

pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

# --. Interference and contact pressure analysis --..
def plot_interference_pressure():
    """Calculate interference and plot contact pressure as function of interference.
    
    Following the formula from coaxial_contact verification:
    gap = g0 + surf(ur, r, bci) - surf(ur, r, b)
    where b is pellet outer radius, bci is clad inner radius
    """
    # Read geometry parameters
    with open(os.path.join(CASE_DIR, "geometry.yaml")) as f:
        geom = yaml.safe_load(f)
    b = float(geom["outer_radius_1"])     # pellet outer radius
    bci = float(geom["inner_radius_2"])  # clad inner radius
    
    # Read contact parameters
    with open(os.path.join(CASE_DIR, "input.yaml")) as f:
        inp = yaml.safe_load(f)
    cc = inp["models"]["contact"]
    g0 = float(cc["initial_gap"])
    k_pen = float(cc["penalty_stiffness"])
    
    # Helper functions for surface average and area-weighted mean
    def surf(v, r, r0):
        return v[np.abs(r - r0) < 2e-5].mean()
    
    def amean(v, r, lo, hi):
        mask = (r >= lo) & (r <= hi)
        return np.sum(v[mask] * r[mask]) / np.sum(r[mask])
    
    # Load history data from CSV
    with open(HISTORY) as f:
        rows = list(csv.DictReader(f))
    
    # Extract data arrays
    time_days = np.array([float(r["time_days"]) for r in rows])
    contact_pressure = np.array([float(r["contact_pressure_MPa"]) for r in rows])
    gap_um = np.array([float(r["gap_um"]) for r in rows])
    
    # Calculate interference from gap data (negative gap = interference)
    # This is equivalent to: gap = g0 + surf(ur, r, bci) - surf(ur, r, b)
    # where negative gap indicates interference
    interference = np.maximum(0, -gap_um * 1e-6)  # convert um to m, take positive part
    
    # Calculate pellet temperature from history data (T_max is used as approximation)
    # In a full implementation, this would use: Tp = amean(T, r, 0.0, b) from field data
    T_pellet = np.array([float(r["T_max_K"]) for r in rows])
    
    # Create plot
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
    
    # Left plot: Contact pressure vs time
    ax1.plot(time_days, contact_pressure, "r-o", lw=2, ms=4, label="Contact pressure")
    ax1.set_xlabel("Time (days)")
    ax1.set_ylabel("Contact pressure (MPa)")
    ax1.set_title("Contact pressure evolution")
    ax1.grid(True, ls=":", alpha=0.6)
    ax1.legend()
    
    # Middle plot: Contact pressure vs interference
    ax2.plot(interference * 1e6, contact_pressure, "b-s", lw=2, ms=4, label="p vs interference")
    ax2.set_xlabel("Interference (μm)")
    ax2.set_ylabel("Contact pressure (MPa)")
    ax2.set_title("Contact pressure as function of interference")
    ax2.grid(True, ls=":", alpha=0.6)
    ax2.legend()
    
    # Right plot: Interference vs pellet temperature
    ax3.plot(T_pellet, interference * 1e6, "g-^", lw=2, ms=4, label="interference vs T")
    ax3.set_xlabel("Pellet temperature (K)")
    ax3.set_ylabel("Interference (μm)")
    ax3.set_title("Interference as function of pellet temperature")
    ax3.grid(True, ls=":", alpha=0.6)
    ax3.legend()
    
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "interference_pressure.png"), dpi=150)
    plt.close(fig)
    print("[INFO] interference_pressure.png saved")
    
    # Print summary statistics
    contact_steps = contact_pressure > 0
    if contact_steps.any():
        print(f"[INFO] Contact initiated at time: {time_days[contact_steps][0]:.1f} days")
        print(f"[INFO] Contact initiated at interference: {interference[contact_steps][0]*1e6:.3f} μm")
        print(f"[INFO] Contact initiated at pellet temperature: {T_pellet[contact_steps][0]:.1f} K")
        print(f"[INFO] Final interference: {interference[-1]*1e6:.3f} μm")
        print(f"[INFO] Final contact pressure: {contact_pressure[-1]:.2f} MPa")
        print(f"[INFO] Final pellet temperature: {T_pellet[-1]:.1f} K")
    else:
        print("[INFO] No contact detected during simulation")


def print_outer_cylinder_mean_stress():
    """Extract the 3x3 stress tensor for each XDMF time step and print its mean."""
    h5_path = XDMF_FILE.replace(".xdmf", ".h5")
    with h5py.File(h5_path, "r") as f:
        if "Function" not in f:
            raise RuntimeError(f"No 'Function' group found in {h5_path}")

        stress_candidates = [
            name for name in f["Function"].keys()
            if "stress" in name.lower()
        ]
        if not stress_candidates:
            raise RuntimeError(
                f"No stress field found in {h5_path}. Available: {list(f['Function'].keys())}"
            )

        stress_field = stress_candidates[0]
        print(f"[INFO] Using stress field '{stress_field}' from {h5_path}")

        def parse_step_name(value):
            try:
                return (0, float(value))
            except Exception:
                return (1, value)

        step_names = sorted(
            list(f["Function"][stress_field].keys()),
            key=parse_step_name,
        )

    r2_i = float(geom["inner_radius_2"])
    r2_o = float(geom["outer_radius_2"])

    for step_index, step_name in enumerate(step_names):
        x_s, y_s, z_s, stress_data = extract_field_xdmf(
            XDMF_FILE, field_name=stress_field, step_index=step_index
        )

        if stress_data.ndim != 2 or stress_data.shape[1] != 9:
            raise RuntimeError(
                f"Expected stress tensor field with 9 components at step {step_name}, got shape {stress_data.shape}"
            )

        r = x_s
        mask = (r >= r2_i - 1e-12) & (r <= r2_o + 1e-12)
        if not np.any(mask):
            print(
                f"[WARNING] No points found in outer cylinder radius range [{r2_i}, {r2_o}] at step {step_name}"
            )
            continue

        stress_tensors = stress_data.reshape(-1, 3, 3)
        outer_stress = stress_tensors[mask]
        mean_stress = outer_stress.mean(axis=0)

        print(f"[INFO] Step {step_name} — Mean stress tensor in outer cylinder (Pa):")
        print(
            f"[[{mean_stress[0,0]: .6e}, {mean_stress[0,1]: .6e}, {mean_stress[0,2]: .6e}]"
            f" [{mean_stress[1,0]: .6e}, {mean_stress[1,1]: .6e}, {mean_stress[1,2]: .6e}]"
            f" [{mean_stress[2,0]: .6e}, {mean_stress[2,1]: .6e}, {mean_stress[2,2]: .6e}]"
        )

plot_interference_pressure()
print_outer_cylinder_mean_stress()

print("\n[INFO] non-regression completed.\n")
