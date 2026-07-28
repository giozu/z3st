#!/usr/bin/env python3
"""
Z3ST study: find the critical internal cracking pressure for the lenticular
cavity case (Neumann traction ramped on region 'cavity').

Strategy: a single ramp from 0 to P_MAX already samples every intermediate
pressure, so instead of re-running many independent simulations to
different pressure caps, we run one coarse, wide-range, few-step ramp to
locate the pressure at which damage stops being smooth/diffuse (D_max
accelerates sharply or the staggered solver stops converging in a couple of
iterations), then a second, fine-step ramp restricted to a narrow window
around that bracket to pin the critical pressure down precisely.
"""

import os
import re
import subprocess
import numpy as np
import meshio
import yaml
import matplotlib.pyplot as plt

CASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(CASE_DIR, "output")
BC_FILE = os.path.join(CASE_DIR, "boundary_conditions.yaml")
INPUT_FILE = os.path.join(CASE_DIR, "input.yaml")


def write_boundary_conditions(p_max, n_steps):
    steps = np.linspace(0, p_max, n_steps).tolist()
    config = {
        'mechanical': {
            'uo2': [
                {'type': 'Clamp_y', 'region': 'ymin'},
                {'type': 'Clamp_x', 'region': 'xmin'},
                {
                    'type': 'Neumann',
                    'region': 'cavity',
                    'traction': steps,
                },
            ]
        }
    }
    with open(BC_FILE, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)


def set_n_steps(n_steps, p_max):
    with open(INPUT_FILE, 'r') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if line.strip().startswith('n_steps:'):
            lines[i] = f"n_steps: {n_steps}  #For {n_steps} steps (0 to {p_max:.3e} Pa)\n"
    with open(INPUT_FILE, 'w') as f:
        f.writelines(lines)


def set_stag_tol(tol):
    """Tighten stag_tol under both the mechanical: and damage: sections.

    Track the current top-level section while scanning so the two
    (differently-scoped) 'stag_tol:' lines are both updated without
    touching thermal's.
    """
    with open(INPUT_FILE, 'r') as f:
        lines = f.readlines()
    section = None
    for i, line in enumerate(lines):
        stripped = line.strip()
        if not line.startswith((' ', '\t')) and stripped.endswith(':'):
            section = stripped[:-1]
        if section in ("mechanical", "damage") and stripped.startswith('stag_tol:'):
            lines[i] = f"  stag_tol: {tol:.1e}\n"
    with open(INPUT_FILE, 'w') as f:
        f.writelines(lines)


def run_sim(log_name):
    if os.path.exists(OUTPUT_DIR):
        for f in os.listdir(OUTPUT_DIR):
            if f.startswith("simulation_") and f.endswith(".vtu"):
                os.remove(os.path.join(OUTPUT_DIR, f))
    else:
        os.makedirs(OUTPUT_DIR)

    env = os.environ.copy()
    z3st_root = os.path.abspath(os.path.join(CASE_DIR, "../../.."))
    env["PYTHONPATH"] = z3st_root + ":" + env.get("PYTHONPATH", "")

    log_path = os.path.join(CASE_DIR, log_name)
    with open(log_path, "w") as log_file:
        subprocess.run(["python3", "-m", "z3st"], cwd=CASE_DIR, check=True,
                        stdout=log_file, stderr=subprocess.STDOUT, env=env)
    return log_path


def parse_iterations(log_path):
    """Return list of (step_index, n_iterations) parsed from the log."""
    with open(log_path) as f:
        text = f.read()
    steps = [int(m) for m in re.findall(r'## Step (\d+)/\d+', text)]
    iters = [int(m) for m in re.findall(r'converged in (\d+) iterations', text)]
    return list(zip(steps, iters))


def extract_damage_curve(p_max, n_steps):
    """Return arrays (pressure, D_max) for every solved step, reading vtu outputs."""
    def step_of(fname):
        m = re.search(r'_(\d+)\.vtu', fname)
        return int(m.group(1)) if m else -1

    vtu_files = sorted(
        (f for f in os.listdir(OUTPUT_DIR) if f.startswith("simulation_") and f.endswith(".vtu")),
        key=step_of,
    )
    pressures, d_max = [], []
    steps_axis = np.linspace(0, p_max, n_steps)
    for fname in vtu_files:
        step = step_of(fname)
        m = meshio.read(os.path.join(OUTPUT_DIR, fname))
        d_max.append(float(np.max(m.point_data["Damage"])))
        pressures.append(steps_axis[step] if 0 <= step < n_steps else np.nan)
    return np.array(pressures), np.array(d_max)


def find_bracket(pressures, d_max, threshold=0.5):
    """First index where D_max crosses `threshold`; bracket is the pressure
    just before and just after that crossing."""
    above = np.where(d_max >= threshold)[0]
    if len(above) == 0:
        return None
    i = above[0]
    if i == 0:
        return (0.0, pressures[0])
    return (pressures[i - 1], pressures[i])


def plot_results(coarse, refine, ultra, p_critical, out_path):
    """coarse/refine/ultra are each (pressures_Pa, D_max) tuples or None if that
    stage never ran (e.g. an earlier stage never crossed the threshold)."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    ax = axes[0]
    if coarse is not None:
        p_c, d_c = coarse
        ax.plot(p_c / 1e6, d_c, 'o-', color='0.4', markersize=3, label="Coarse sweep")
    if refine is not None:
        p_r, d_r = refine
        ax.axvspan(p_r[0] / 1e6, p_r[-1] / 1e6, color='tab:orange', alpha=0.15,
                   label="Refine window")
    ax.axhline(0.5, color='gray', ls='--', lw=1, alpha=0.7)
    ax.set_xlabel("Applied cavity pressure (MPa)")
    ax.set_ylabel("max(Damage) (-)")
    ax.set_title("Coarse sweep (full range)")
    ax.grid(True, ls=':', alpha=0.5)
    ax.legend(fontsize=8)

    ax = axes[1]
    if refine is not None:
        p_r, d_r = refine
        ax.plot(p_r / 1e6, d_r, 'o-', color='tab:orange', markersize=3, label="Refine sweep")
    if ultra is not None:
        p_u, d_u = ultra
        ax.plot(p_u / 1e6, d_u, 'o-', color='tab:red', markersize=3, label="Ultra-fine sweep")
    ax.axhline(0.5, color='gray', ls='--', lw=1, alpha=0.7, label="D = 0.5 threshold")
    if p_critical is not None:
        ax.axvline(p_critical / 1e6, color='k', ls=':', lw=1.5,
                   label=f"Critical pressure ~ {p_critical / 1e6:.2f} MPa")
    ax.set_xlabel("Applied cavity pressure (MPa)")
    ax.set_ylabel("max(Damage) (-)")
    ax.set_title("Refine / ultra-fine sweeps (zoomed on transition)")
    ax.grid(True, ls=':', alpha=0.5)
    ax.legend(fontsize=8)

    fig.suptitle("Cavity pressurization: damage evolution across sweep stages")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"[INFO] Plot saved: {out_path}")


def main():
    print("=" * 60)
    print(" Coarse pressure sweep: locating the critical bracket")
    print("=" * 60)
    P_MAX_COARSE = 1.5e9   # 1.5 GPa, 10x the placeholder guess
    N_STEPS_COARSE = 150

    write_boundary_conditions(P_MAX_COARSE, N_STEPS_COARSE)
    set_n_steps(N_STEPS_COARSE, P_MAX_COARSE)
    log_path = run_sim("log_study_coarse.txt")

    pressures, d_max = extract_damage_curve(P_MAX_COARSE, N_STEPS_COARSE)
    iters = parse_iterations(log_path)
    max_iters_used = max(n for _, n in iters) if iters else None

    print(f"Coarse sweep done. Max iterations used at any step: {max_iters_used}")
    print(f"D_max final: {d_max[-1]:.4f} at P = {pressures[-1]:.3e} Pa")

    bracket = find_bracket(pressures, d_max, threshold=0.5)
    if bracket is None:
        print("\n[WARNING] D_max never reached 0.5 within the coarse range.")
        print("Increase P_MAX_COARSE and re-run before attempting the refine stage.")
        np.savetxt(os.path.join(CASE_DIR, "pressure_sweep_coarse.csv"),
                   np.column_stack([pressures, d_max]),
                   header="pressure_Pa,D_max", delimiter=",", comments="")
        plot_results((pressures, d_max), None, None, None,
                     os.path.join(CASE_DIR, "pressure_sweep_summary.png"))
        return

    print(f"\nBracket for critical pressure: {bracket[0]:.3e} Pa -- {bracket[1]:.3e} Pa")
    np.savetxt(os.path.join(CASE_DIR, "pressure_sweep_coarse.csv"),
               np.column_stack([pressures, d_max]),
               header="pressure_Pa,D_max", delimiter=",", comments="")

    print("\n" + "=" * 60)
    print(" Refine sweep: fine ramp restricted to the bracket")
    print("=" * 60)
    P_MAX_REFINE = bracket[1] * 1.1  # small margin past the upper bracket
    N_STEPS_REFINE = 401

    write_boundary_conditions(P_MAX_REFINE, N_STEPS_REFINE)
    set_n_steps(N_STEPS_REFINE, P_MAX_REFINE)
    log_path = run_sim("log_study_refine.txt")

    pressures_r, d_max_r = extract_damage_curve(P_MAX_REFINE, N_STEPS_REFINE)
    np.savetxt(os.path.join(CASE_DIR, "pressure_sweep_refine.csv"),
               np.column_stack([pressures_r, d_max_r]),
               header="pressure_Pa,D_max", delimiter=",", comments="")

    iters_r = parse_iterations(log_path)
    max_iters_refine = max(n for _, n in iters_r) if iters_r else None
    print(f"Refine sweep done. Max iterations used at any step: {max_iters_refine}")

    bracket_r = find_bracket(pressures_r, d_max_r, threshold=0.5)
    if not bracket_r:
        print("\n[WARNING] Refine ramp did not cross D_max=0.5 either -- widen P_MAX_REFINE.")
        plot_results((pressures, d_max), (pressures_r, d_max_r), None, None,
                     os.path.join(CASE_DIR, "pressure_sweep_summary.png"))
        return

    p_c_refine = 0.5 * (bracket_r[0] + bracket_r[1])
    print(f"\n[RESULT] Refine estimate (D_max=0.5 crossing): ~{p_c_refine:.4e} Pa")
    print(f"Bracket for ultra-fine pass: {bracket_r[0]:.4e} Pa -- {bracket_r[1]:.4e} Pa")

    print("\n" + "=" * 60)
    print(" Ultra-fine sweep: zoomed on the refine bracket, tightened stag_tol")
    print("=" * 60)
    # Same numerical-artifact risk as coarse->refine can in principle recur
    # if the refine step size was still too coarse right at the transition,
    # so this pass both zooms in further AND tightens the staggering
    # tolerance (default 1e-4 -> 1e-5) to make the alternating solve track
    # the true equilibrium path more faithfully, rather than declaring
    # convergence on a looser "close enough between iterations" criterion.
    STAG_TOL_ULTRA = 1e-5
    P_MAX_ULTRA = bracket_r[1] * 1.05  # small margin past the refined upper bracket
    N_STEPS_ULTRA = 401

    write_boundary_conditions(P_MAX_ULTRA, N_STEPS_ULTRA)
    set_n_steps(N_STEPS_ULTRA, P_MAX_ULTRA)
    set_stag_tol(STAG_TOL_ULTRA)
    log_path = run_sim("log_study_ultra.txt")

    pressures_u, d_max_u = extract_damage_curve(P_MAX_ULTRA, N_STEPS_ULTRA)
    np.savetxt(os.path.join(CASE_DIR, "pressure_sweep_ultra.csv"),
               np.column_stack([pressures_u, d_max_u]),
               header="pressure_Pa,D_max", delimiter=",", comments="")

    iters_u = parse_iterations(log_path)
    max_iters_ultra = max(n for _, n in iters_u) if iters_u else None
    print(f"Ultra-fine sweep done. Max iterations used at any step: {max_iters_ultra}")
    print("(compare to the coarse sweep's max iterations above -- a much lower "
          "value here supports that this pass is tracking the true equilibrium "
          "path rather than struggling near a spurious staggered minimum.)")

    bracket_u = find_bracket(pressures_u, d_max_u, threshold=0.5)
    p_c = None
    if bracket_u:
        p_c = 0.5 * (bracket_u[0] + bracket_u[1])
        print(f"\n[RESULT] Critical cracking pressure (D_max=0.5 crossing): ~{p_c:.4e} Pa "
              f"(bracket {bracket_u[0]:.4e} -- {bracket_u[1]:.4e} Pa)")
    else:
        print("\n[WARNING] Ultra-fine ramp did not cross D_max=0.5 -- widen its P_MAX margin.")

    plot_results((pressures, d_max), (pressures_r, d_max_r), (pressures_u, d_max_u), p_c,
                 os.path.join(CASE_DIR, "pressure_sweep_summary.png"))


if __name__ == "__main__":
    main()
