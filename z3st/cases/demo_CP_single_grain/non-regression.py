#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: demo_CP_single_grain

non-regression script
---------------------
Single crystal plasticity with automatic differentiation.
Uniaxial tension test on a single grain cube.

"""

import os
import glob
import yaml
import numpy as np
import matplotlib.pyplot as plt

from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_verification import *

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
OUTPUT_DIR = os.path.join(CASE_DIR, "output")
OUT_JSON = os.path.join(OUTPUT_DIR, "non-regression.json")
INPUT_FILE = os.path.join(CASE_DIR, "input.yaml")

# Load input configuration
with open(INPUT_FILE, 'r') as f:
    input_data = yaml.safe_load(f)

# Load material data
mat_path = os.path.join(CASE_DIR, input_data['materials']['grain'])
with open(mat_path, 'r') as f:
    mat_data = yaml.safe_load(f)

E = float(mat_data.get('E', 200e9))
nu = float(mat_data.get('nu', 0.3))
g0 = float(mat_data.get('g0', 50e6))
gamma0 = float(mat_data.get('gamma0', 0.02))
n_pow = float(mat_data.get('n_pow', 10.0))

print(f"[INFO] Material loaded: E = {E:.2e} Pa, nu = {nu}")
print(f"[INFO] Plasticity params: g0 = {g0*1e-6:.1f} MPa, gamma0 = {gamma0}, n = {n_pow}")

APPLIED_STRAIN_MAX = 0.01  # from boundary conditions (10 mm displacement on 1 m cube)
TOLERANCE = 0.25  # relaxed tolerance for plasticity (25%)

# --.. ..- .-.. .-.. --- analytic functions --.. ..- .-.. .-.. ---
def elastic_stress(epsilon, E, nu):
    """Analytical elastic stress for constrained uniaxial loading (oedometric)."""
    M = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu))
    return M * epsilon


def viscoplastic_saturation_stress(strain_rate, schmid_factor, g0, gamma0, n_pow):
    """
    Analytical saturation stress for viscoplastic power law.

    For a viscoplastic material without hardening under constant strain rate,
    the stress saturates when plastic strain rate equals total strain rate.

    At steady state (dσ/dt = 0):
        ε̇_total = ε̇_plastic
        ε̇_total = m · γ̇ = m · γ₀ · (m·σ/g₀)ⁿ

    Solving for σ:
        σ_sat = (g₀/m) · (ε̇_total / (m·γ₀))^(1/n)

    Parameters:
        strain_rate: Total strain rate (1/s)
        schmid_factor: Schmid factor m (dimensionless)
        g0: Slip resistance (Pa)
        gamma0: Reference slip rate (1/s)
        n_pow: Power law exponent

    Returns:
        σ_sat: Saturation stress (Pa)
    """
    m = schmid_factor
    sigma_sat = (g0 / m) * (strain_rate / (m * gamma0))**(1.0 / n_pow)
    return sigma_sat


# --.. ..- .-.. .-.. --- checks --.. ..- .-.. .-.. ---
vtu_files = sorted(glob.glob(os.path.join(OUTPUT_DIR, "fields_*.vtu")))
print(f"[INFO] Processing {len(vtu_files)} VTU files...")

if len(vtu_files) == 0:
    raise FileNotFoundError(f"No VTU files found in {OUTPUT_DIR}")

# List available fields in last VTU
print(f"[INFO] Available fields in {os.path.basename(vtu_files[-1])}:")
list_fields(vtu_files[-1])

# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
strains_np = []
stresses_np = []

for vtu in vtu_files:
    # Extract displacement field
    _, _, _, disp = extract_field(vtu, "Displacement")

    # Calculate strain from displacement (geometry: 1×1×1 cube)
    uz_max = np.max(disp[:, 2])
    uz_min = np.min(disp[:, 2])
    epsilon_zz = (uz_max - uz_min) / 1.0  # Lz = 1.0 m
    strains_np.append(epsilon_zz)

    # Extract stress field (cell data)
    _, _, _, stress = extract_field(vtu, "Stress_grain (cells)")
    sigma_zz = np.mean(stress[:, 8])
    stresses_np.append(sigma_zz)

strains_np = np.array(strains_np)
stresses_np = np.array(stresses_np)

# Analytical elastic prediction
elastic_stresses_np = elastic_stress(strains_np, E, nu)

# --.. ..- .-.. .-.. --- semi-analytical saturation stress --.. ..- .-.. .-.. ---
# Calculate strain rate (assuming constant rate loading)
time_total = float(input_data['time'][-1])  # Final time
strain_rate = APPLIED_STRAIN_MAX / time_total  # ε̇ = Δε/Δt

# Schmid factor for (111)[0-11] slip system under z-loading
schmid_factor = 0.408248  # From Schmid factor calculation

# Calculate analytical saturation stress
sigma_sat_analytical = viscoplastic_saturation_stress(
    strain_rate=strain_rate,
    schmid_factor=schmid_factor,
    g0=g0,
    gamma0=gamma0,
    n_pow=n_pow
)

print(f"\n[ANALYTICAL SOLUTION]")
print(f"  Strain rate ε̇_total:     {strain_rate:.4e} s⁻¹")
print(f"  Schmid factor m:         {schmid_factor:.6f}")
print(f"  Saturation stress σ_sat: {sigma_sat_analytical*1e-6:.2f} MPa")
print(f"  (At steady state: ε̇_plastic = ε̇_total)")

# --.. ..- .-.. .-.. --- plotting --.. ..- .-.. .-.. ---
Pa_to_MPa = 1e-6

fig, ax = plt.subplots(figsize=(10, 7))

# Plot stress-strain curves
ax.plot(strains_np * 100, stresses_np * Pa_to_MPa, 'o-',
        label='Z3ST (Crystal Plasticity)', color='#E63946',
        linewidth=2.5, markersize=6, alpha=0.8)
ax.plot(strains_np * 100, elastic_stresses_np * Pa_to_MPa, '--',
        label='Elastic Reference', color='#1D3557',
        linewidth=2, alpha=0.7)

# Analytical saturation stress (horizontal asymptote)
ax.axhline(y=sigma_sat_analytical * Pa_to_MPa, color='#2A9D8F', linestyle='-.',
           linewidth=2.5, label=f'Analytical Saturation σ_sat = {sigma_sat_analytical*Pa_to_MPa:.1f} MPa',
           alpha=0.8)

# Estimate yield stress (approximate, based on Schmid factor)
sigma_yield_approx = g0 / schmid_factor
ax.axhline(y=sigma_yield_approx * Pa_to_MPa, color='gray', linestyle=':',
           linewidth=1.5, label=f'Approx. Yield ({sigma_yield_approx*Pa_to_MPa:.0f} MPa)')

ax.set_xlabel(r'Strain $\varepsilon_{zz}$ [%]', fontsize=13, fontweight='bold')
ax.set_ylabel(r'Stress $\sigma_{zz}$ [MPa]', fontsize=13, fontweight='bold')
ax.set_title('Crystal Plasticity: Uniaxial Loading (Single Grain)',
             fontsize=14, fontweight='bold', pad=15)
ax.grid(True, alpha=0.3, linestyle='--')
ax.legend(fontsize=11, frameon=True, loc='best')

plt.tight_layout()
plot_path = os.path.join(OUTPUT_DIR, "stress_strain_curve.png")
plt.savefig(plot_path, dpi=300)
print(f"[INFO] Plot saved in: {plot_path}")

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---

# Final state values (last time step)
epsilon_zz_final = strains_np[-1]
sigma_zz_final = stresses_np[-1]
sigma_elastic_final = elastic_stresses_np[-1]

# Calculate stress reduction (should be negative for plasticity)
stress_reduction = (sigma_zz_final - sigma_elastic_final) / sigma_elastic_final
rel_error_sigma = abs(stress_reduction)

# Check strain reached target
rel_error_strain = abs(epsilon_zz_final - APPLIED_STRAIN_MAX) / APPLIED_STRAIN_MAX

# Check convergence to saturation stress
saturation_error = abs(sigma_zz_final - sigma_sat_analytical) / sigma_sat_analytical

print(f"\n[RESULTS] Final State (Step {len(vtu_files)-1}):")
print(f"  Strain ε_zz:        {epsilon_zz_final:.4e} (target: {APPLIED_STRAIN_MAX:.4e})")
print(f"  Stress σ_zz (CP):   {sigma_zz_final*Pa_to_MPa:.2f} MPa")
print(f"  Elastic prediction: {sigma_elastic_final*Pa_to_MPa:.2f} MPa")
print(f"  Analytical σ_sat:   {sigma_sat_analytical*Pa_to_MPa:.2f} MPa")
print(f"  Stress reduction:   {stress_reduction*100:+.1f}%")
print(f"  Saturation error:   {saturation_error*100:+.1f}%")

# Check if stress is approaching saturation
print(f"\n[VERIFICATION] Semi-Analytical Check:")
if saturation_error < 0.10:  # Within 10%
    print(f"  ✓ EXCELLENT: Stress converged to analytical saturation (error = {saturation_error*100:.1f}%)")
elif saturation_error < 0.25:  # Within 25%
    print(f"  ✓ GOOD: Stress approaching saturation (error = {saturation_error*100:.1f}%)")
else:
    print(f"  ✗ WARNING: Stress deviates from saturation (error = {saturation_error*100:.1f}%)")
    print(f"    This may indicate: insufficient loading, wrong parameters, or numerical issues.")

# Error dictionary
errors = {
    "sigma_zz_final": {
        "numerical": float(sigma_zz_final),
        "reference": float(sigma_sat_analytical),  # Compare vs σ_sat, not elastic!
        "abs_error": float(abs(sigma_zz_final - sigma_sat_analytical)),
        "rel_error": float(saturation_error),  # Already calculated above
    },
    "epsilon_zz_final": {
        "numerical": float(epsilon_zz_final),
        "reference": float(APPLIED_STRAIN_MAX),
        "abs_error": float(abs(epsilon_zz_final - APPLIED_STRAIN_MAX)),
        "rel_error": float(rel_error_strain),
    },
    "saturation_convergence": {
        "numerical": float(saturation_error * 100),  # Error percentage
        "reference": 0.0,  # Perfect convergence target
        "abs_error": float(saturation_error * 100),
        "rel_error": float(saturation_error),  # Relative to perfect convergence
    },
}

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] non-regression completed.\n")
