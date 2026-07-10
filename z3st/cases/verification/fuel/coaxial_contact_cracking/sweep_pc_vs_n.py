#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST analysis script --.. ..- .-.. .-.. ---
"""
Z3ST case: verification/fuel/coaxial_contact_cracking

Sweeps the rod-average LHR over a fixed grid, running a fresh, independent
z3st case for every value (each ramping 0 -> LHR_target over N_STEPS_PER_RUN
steps), and records the converged contact pressure Pc at each run's own LAST
time step -- once with the isotropic-softening cracking model active
(models/cracking_model.py) and once with it switched off, so the softening
effect can be read directly off one figure.

The two variants share every input file except fuel.yaml's "cracking" key:
the "cracking" run uses this case's own fuel.yaml (cracking: isotropic), the
"no cracking" run uses the same material with that key dropped -- both are
generated in-memory, so this case has no dependency on coaxial_contact.

n (number of macro-cracks) is not an exported field, so it is parsed from
the last "[cracking]" line of each cracking run's log_z3st.md (the converged
sample of that run's last step). Pc is read directly from the
"ContactPressure" cell field of each run's last VTU output. The no-cracking
curve is plotted against the CRACKING run's own n(LHR) -- n has no meaning
for the no-cracking material (it stays virgin), so this is used purely as a
shared x-axis tied to LHR, letting the two curves be compared at matching
operating points.
"""

import os
import re
import csv
import glob
import shutil
import subprocess
import sys

import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import yaml

CASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(CASE, "output")
SWEEP_DIR = os.path.join(CASE, "sweep_runs")

PYTHON = sys.executable
SHARED_FILES = ["geometry.yaml", "boundary_conditions.yaml", "clad.yaml", "mesh.msh"]

LHR_MIN, LHR_MAX, LHR_STEP = 0.0, 100_000.0, 10_000.0
N_STEPS_PER_RUN = 13   # ramp resolution within each independent run

CRACK_RE = re.compile(
    r"\[cracking\]\s+\w+:\s+LHR_max\s+=\s+([\d.]+)\s+kW/m.*?"
    r"n\s+=\s+([\d.]+)\s+cracks,\s+E_iso/E\s+=\s+([\d.]+)"
)

base_input = yaml.safe_load(open(os.path.join(CASE, "input.yaml")))
base_fuel = yaml.safe_load(open(os.path.join(CASE, "fuel.yaml")))

fuel_variants = {
    "cracking": dict(base_fuel),
    "nocracking": {k: v for k, v in base_fuel.items() if k != "cracking"},
}


def run_at(lhr_target, variant):
    """Run a fresh, independent case (fuel material = fuel_variants[variant])
    ramping 0 -> lhr_target and return (n, E_iso_over_E, Pc_MPa) at its own
    last time step. n/E_iso_over_E are 0/1 (virgin) for the no-cracking
    variant, since the cracking model never activates."""
    rundir = os.path.join(SWEEP_DIR, f"lhr_{int(lhr_target):07d}_{variant}")
    os.makedirs(rundir, exist_ok=True)
    for fname in SHARED_FILES:
        shutil.copy(os.path.join(CASE, fname), os.path.join(rundir, fname))
    with open(os.path.join(rundir, "fuel.yaml"), "w") as f:
        yaml.safe_dump(fuel_variants[variant], f, sort_keys=False)

    run_input = dict(base_input)
    run_input["lhr"] = [0.0, float(lhr_target)]
    run_input["time"] = [0.0, 1.0]
    run_input["n_steps"] = N_STEPS_PER_RUN
    with open(os.path.join(rundir, "input.yaml"), "w") as f:
        yaml.safe_dump(run_input, f, sort_keys=False)

    log_path = os.path.join(rundir, "log_z3st.md")
    with open(log_path, "w") as log_f:
        subprocess.run(
            [PYTHON, "-m", "z3st"], cwd=rundir,
            stdout=log_f, stderr=subprocess.STDOUT, check=True,
        )

    vtu_files = sorted(glob.glob(os.path.join(rundir, "output", "fields_*.vtu")))
    if not vtu_files:
        raise RuntimeError(f"no VTU output found for LHR={lhr_target} ({variant}, see {log_path})")
    m = pv.read(vtu_files[-1])
    pc = float(m.cell_data["ContactPressure"][0]) / 1e6   # MPa

    if variant == "nocracking":
        return 0.0, 1.0, pc

    log_text = open(log_path).read()
    matches = CRACK_RE.findall(log_text)
    if not matches:
        raise RuntimeError(f"no '[cracking]' line found for LHR={lhr_target} (see {log_path})")
    _, n, er = matches[-1]   # converged (last) sample of the run's last step
    return float(n), float(er), pc


lhr_grid = np.arange(LHR_MIN, LHR_MAX + LHR_STEP / 2, LHR_STEP)
n_cracks, E_ratio, p_crack, p_nocrack = [], [], [], []

for lhr in lhr_grid:
    print(f"[INFO] LHR = {lhr/1e3:.1f} kW/m ...")
    n, er, pc_c = run_at(lhr, "cracking")
    _, _, pc_nc = run_at(lhr, "nocracking")
    n_cracks.append(n)
    E_ratio.append(er)
    p_crack.append(pc_c)
    p_nocrack.append(pc_nc)
    print(f"       n = {n:.4f}  E_iso/E = {er:.6f}  Pc(cracking) = {pc_c:.4f} MPa  Pc(no cracking) = {pc_nc:.4f} MPa")

lhr_grid, n_cracks, E_ratio, p_crack, p_nocrack = map(
    np.array, (lhr_grid, n_cracks, E_ratio, p_crack, p_nocrack))

# --. table --..
print(f"\n{'LHR (kW/m)':>10}  {'n':>8}  {'E_iso/E':>9}  {'Pc crack (MPa)':>15}  {'Pc no-crack (MPa)':>18}")
for lhr, n, er, pc_c, pc_nc in zip(lhr_grid, n_cracks, E_ratio, p_crack, p_nocrack):
    print(f"{lhr/1e3:10.1f}  {n:8.4f}  {er:9.6f}  {pc_c:15.4f}  {pc_nc:18.4f}")

csv_path = os.path.join(OUT, "sweep_pc_vs_n.csv")
os.makedirs(OUT, exist_ok=True)
with open(csv_path, "w", newline="") as fh:
    writer = csv.writer(fh)
    writer.writerow(["lhr_W_per_m", "n_cracks", "E_iso_over_E", "Pc_cracking_MPa", "Pc_nocracking_MPa"])
    for lhr, n, er, pc_c, pc_nc in zip(lhr_grid, n_cracks, E_ratio, p_crack, p_nocrack):
        writer.writerow([lhr, n, er, pc_c, pc_nc])
print(f"\n[INFO] {csv_path} written")

# --. Pc vs n plot (both curves share the cracking run's n(LHR) x-axis) --..
fig, ax1 = plt.subplots(figsize=(8, 6))
ax1.plot(n_cracks, p_crack, "o-", color="C3", lw=2, label="Pc, cracking")
ax1.plot(n_cracks, p_nocrack, "s--", color="C0", lw=2, label="Pc, no cracking")
ax1.set_xlabel("number of cracks n (-)")
ax1.set_ylabel("contact pressure Pc (MPa)")
ax1.set_title("Contact pressure vs. number of cracks (LHR sweep, 0-100 kW/m)")
ax1.grid(True, ls=":", alpha=0.6)

# Every LHR value is still labelled, but points bunch up tightly in n as
# cracking saturates, so text placed right next to each point overlaps. Fan
# the labels out at evenly-spaced heights (in axes-fraction y, ordered by n)
# with a thin leader line back to the real data point instead.
order = np.argsort(n_cracks)
y_fracs = np.linspace(0.06, 0.94, len(order))
for rank, idx in enumerate(order):
    x = n_cracks[idx]
    y_anchor = max(p_crack[idx], p_nocrack[idx])
    ax1.annotate(
        f"{lhr_grid[idx] / 1e3:.0f} kW/m",
        xy=(x, y_anchor), xycoords="data",
        xytext=(x, y_fracs[rank]), textcoords=("data", "axes fraction"),
        fontsize=7.5, color="0.3", ha="center",
        arrowprops=dict(arrowstyle="-", lw=0.6, color="0.7", shrinkA=0, shrinkB=2),
    )
ax1.legend(loc="upper left")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "sweep_pc_vs_n.png"), dpi=150)
print("[INFO] sweep_pc_vs_n.png saved")

print("\n[INFO] sweep_pc_vs_n analysis completed.\n")
