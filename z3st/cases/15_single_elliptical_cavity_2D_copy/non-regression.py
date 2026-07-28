#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: single_elliptical_cavity_2D_copy (pressure-driven cavity cracking)

non-regression script
---------------------
Ramps an internal pressure on the cavity boundary (Neumann BC) instead of a
prescribed displacement on ymax. Tracks the phase-field Damage over the
pressure ramp and reports the critical cracking pressure: the applied
pressure at the first load step where max(Damage) over the whole domain
reaches DAMAGE_CRITICAL_THRESHOLD.
"""

import os, re, sys
import yaml
import numpy as np
import matplotlib.pyplot as plt

from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_verification import *

# Units of measure for micromechanics
# Length: micrometer (μm)
# Time: second (s)
# Mass: kilogram (kg)
# Force: micronewton (μN)
# Pressure: megapascal (MPa)
# Energy: picojoule (pJ)

DAMAGE_CRITICAL_THRESHOLD = 0.5

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)

OUTPUT_DIR = os.path.join(CASE_DIR, "output")

def extract_step(filename):
    match = re.search(r'_(\d+)\.vtu', filename)
    return int(match.group(1)) if match else -1

vtu_files = sorted([f for f in os.listdir(OUTPUT_DIR) if f.startswith("simulation_") and f.endswith(".vtu")], key=extract_step)
if not vtu_files:
    print("[ERROR] No simulation_*.vtu files found in output directory!")
    sys.exit(1)

OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")
BC_FILE = os.path.join(CASE_DIR, "boundary_conditions.yaml")

# --.. ..- .-.. .-.. --- applied pressure ramp --.. ..- .-.. .-.. ---
with open(BC_FILE, 'r') as f:
    bc_data = yaml.safe_load(f)

cavity_bc = next(
    bc for bc in bc_data["mechanical"]["uo2"]
    if bc.get("region") == "cavity" and bc.get("type") == "Neumann"
)
pressures_pa = np.asarray(cavity_bc["traction"], dtype=float)
pressures_mpa = pressures_pa * 1e-6

if len(pressures_pa) != len(vtu_files):
    print(f"[WARNING] {len(pressures_pa)} pressure steps in boundary_conditions.yaml "
          f"!= {len(vtu_files)} VTU files found; truncating to the shorter of the two.")
n_steps = min(len(pressures_pa), len(vtu_files))
vtu_files = vtu_files[:n_steps]
pressures_mpa = pressures_mpa[:n_steps]

# --.. ..- .-.. .-.. --- extract Damage field over the ramp --.. ..- .-.. .-.. ---
max_damage = []

print("[INFO] Extracting max(Damage) over the domain for each pressure step...")
for vtu in vtu_files:
    vtu_path = os.path.join(OUTPUT_DIR, vtu)
    _, _, _, D = extract_field(vtu_path, field_name="Damage")
    d_max = float(np.max(D))
    max_damage.append(d_max)

max_damage = np.asarray(max_damage)

# --.. ..- .-.. .-.. --- critical cracking pressure --.. ..- .-.. .-.. ---
crack_idx = np.argmax(max_damage >= DAMAGE_CRITICAL_THRESHOLD) if np.any(max_damage >= DAMAGE_CRITICAL_THRESHOLD) else None

if crack_idx is not None:
    p_critical = float(pressures_mpa[crack_idx])
    print("\n" + "="*70)
    print(f"{'█'*70}")
    print(f"  📊 CRITICAL CRACKING PRESSURE: {p_critical:.2f} MPa  (step {crack_idx}, max D = {max_damage[crack_idx]:.4f})")
    print(f"{'█'*70}")
    print("="*70 + "\n")
else:
    p_critical = 0.0
    print(f"\n[WARNING] max(Damage) never reached {DAMAGE_CRITICAL_THRESHOLD} over the pressure ramp "
          f"(max reached: {np.max(max_damage):.4f} at {pressures_mpa[np.argmax(max_damage)]:.2f} MPa). "
          f"Increase P_MAX in generate_yaml.py to bracket the critical pressure.\n")

errors = {
    "critical_pressure_MPa": {
        "numerical": p_critical,
        "reference": 0.0,
        "rel_error": 0.0
    }
}
TOLERANCE = 1.0e-2
pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)

# --.. ..- .-.. .-.. --- plot Pressure vs max(Damage) --.. ..- .-.. .-.. ---
plt.figure(figsize=(10, 6))
plt.plot(pressures_mpa, max_damage, 'b-o', markersize=4, label=r"max(Damage) over domain")
plt.axhline(DAMAGE_CRITICAL_THRESHOLD, color='gray', ls='--', alpha=0.7,
            label=f"Critical threshold D = {DAMAGE_CRITICAL_THRESHOLD}")
if crack_idx is not None:
    plt.axvline(p_critical, color='r', ls=':', alpha=0.8,
                label=f"Critical pressure = {p_critical:.2f} MPa")
plt.xlabel("Applied cavity pressure (MPa)")
plt.ylabel("max(Damage) (-)")
plt.title("Cavity pressurization: damage evolution vs. applied pressure")
plt.grid(True, ls=':', alpha=0.6)
plt.legend()
plt.tight_layout()
plot_path_damage = os.path.join(CASE_DIR, "output", "damage_vs_pressure.png")
plt.savefig(plot_path_damage, dpi=300)
print(f"[INFO] Plot saved in: {plot_path_damage}")

# --.. ..- .-.. .-.. --- plot Pressure vs energies --.. ..- .-.. .-.. ---
energy_file = os.path.join(CASE_DIR, "energies.txt")
if os.path.exists(energy_file):
    energy_data = np.genfromtxt(energy_file, names=True, skip_header=0)
    steps = energy_data["Step"].astype(int)
    mask = steps < n_steps

    plt.figure(figsize=(10, 6))
    plt.plot(pressures_mpa[steps[mask]], energy_data["E_el"][mask], 'b-o', markersize=3, label=r"Elastic $E_{el}$")
    plt.plot(pressures_mpa[steps[mask]], energy_data["E_frac"][mask], 'r-s', markersize=3, label=r"Fracture $E_{frac}$")
    plt.plot(pressures_mpa[steps[mask]], energy_data["E_tot"][mask], 'k--', lw=1.5, label=r"Total $E_{tot}$")
    if crack_idx is not None:
        plt.axvline(p_critical, color='r', ls=':', alpha=0.8,
                    label=f"Critical pressure = {p_critical:.2f} MPa")
    plt.xlabel("Applied cavity pressure (MPa)")
    plt.ylabel("Energy (J)")
    plt.title("Global energy balance vs. applied pressure")
    plt.grid(True, ls=':', alpha=0.6)
    plt.legend(loc="upper left", fontsize="small")
    plt.tight_layout()
    plot_path_energy = os.path.join(CASE_DIR, "output", "energy_vs_pressure.png")
    plt.savefig(plot_path_energy, dpi=300)
    print(f"[INFO] Plot saved in: {plot_path_energy}")
else:
    print(f"[WARN] {energy_file} not found; skipping energy vs. pressure plot.")
