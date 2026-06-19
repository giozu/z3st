#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: plasticity_2D

non-regression script
---------------------

"""

import os
import re
import yaml
import numpy as np
from glob import glob
import matplotlib.pyplot as plt
from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_verification import *

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
OUTPUT_DIR = os.path.join(CASE_DIR, "output")
VTU_FILES = sorted(glob(os.path.join(OUTPUT_DIR, "fields_*.vtu")))
OUT_JSON = os.path.join(OUTPUT_DIR, "non-regression.json")
MATERIAL_FILE = os.path.join(CASE_DIR, "../../../../materials/steel.yaml")
GEOMETRY_FILE = os.path.join(CASE_DIR, "geometry.yaml")
BC_FILE = os.path.join(CASE_DIR, "boundary_conditions.yaml")
MESH_GEO_FILE = os.path.join(CASE_DIR, "mesh.geo")

# Geometry, material, boundary conditions:
with open(GEOMETRY_FILE, 'r') as f:
    geom_data = yaml.safe_load(f)

Lx = float(geom_data.get('Lx'))
Ly = float(geom_data.get('Ly'))

with open(MATERIAL_FILE, 'r') as f:
    mat_data = yaml.safe_load(f)

E     = float(mat_data.get('E'))
nu    = float(mat_data.get('nu'))
k     = float(mat_data.get('k'))
alpha = float(mat_data.get('alpha'))
rho   = float(mat_data.get('rho'))
mu    = float(mat_data.get('mu_gamma'))
q0    = float(mat_data.get('gamma_heating'))
sigma_y = float(mat_data.get('yield_strength'))
H     = float(mat_data.get('hardening_modulus'))

with open(MESH_GEO_FILE, 'r') as f:
    content = f.read()
nxy = int(re.search(r'nxy\s*=\s*(\d+);', content).group(1)) - 1
print(f"[INFO] nFE loaded: nxy = {nxy}")

y_target, mask_tol = (Ly / 2, Ly / (2 * nxy))  # m, m (extraction line selection and tolerance)

# --.. ..- .-.. .-.. --- analytic functions  --.. ..- .-.. .-.. ---
def analytic_stress(epsilon):
    """
    Compute uniaxial stress for a given strain using J2 plasticity with isotropic hardening
    """
    # Effective properties for plane strain uniaxial loading
    E_eff = E / (1 - nu**2)
    phi = np.sqrt(1 - nu + nu**2)
    sigma_y_eff = sigma_y / phi
    eps_y_eff = sigma_y_eff / E_eff
    
    # Elastic regime
    if epsilon <= eps_y_eff:
        return E_eff * epsilon
    
    # Plastic regime (linear isotropic hardening H)
    H_eff = H / (phi**2)
    Et_eff = (E_eff * H_eff) / (E_eff + H_eff)
    
    delta_eps = epsilon - eps_y_eff
    sigma = sigma_y_eff + Et_eff * delta_eps

    return sigma

U_X_REF = np.linspace(0.0, 0.0004, 21)
STRAINS_REF = U_X_REF /Lx
STRESSES_REF_ANALYTIC = [analytic_stress(eps) for eps in STRAINS_REF]

TOLERANCE = 3e-2

# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
print(f"[INFO] Target y-plane for extraction: y = {y_target:.4e} m")

# stress–strain arrays
strains = []
stresses = []
displacements = []
plastic_strains = []

for step, vtufile in enumerate(VTU_FILES):
    print(f"\n[STEP {step}] Processing {os.path.basename(vtufile)}")

    # Stress extraction
    x_S, y_S, z_S, S_all = extract_field(vtufile, field_name="Stress (cells)")
         
    # Check shape/dims
    mask = (np.abs(y_S - y_target) < mask_tol)
    s_val = float(np.mean(S_all[mask, 0]))
    stresses.append(s_val)

    # Displacement extraction
    x_u_all, y_u_all, z_u_all, u_all = extract_field(vtufile, field_name="Displacement")
    mask_u = np.abs(y_u_all - y_target) < mask_tol
    u_max = np.max(u_all[mask_u, 0])
    displacements.append(u_max)

    # Strain extraction
    x_eps_all, y_eps_all, z_eps_all, E_all = extract_field(vtufile, field_name="Strain (cells)")
    mask_e = np.abs(y_eps_all - y_target) < mask_tol
    eps_eng = float(np.mean(E_all[mask_e, 0]))
    strains.append(eps_eng)
    
    # Plastic strain extraction
    x_p, y_p, z_p, P_all = extract_field(vtufile, field_name="CumulativePlasticStrain")
    p_val = float(np.mean(P_all[mask]))
    plastic_strains.append(p_val)

# Curve data
strain_ref_np = np.array(STRAINS_REF, dtype=float)
stresses_analytic_np = np.array(STRESSES_REF_ANALYTIC, dtype=float)
displ_x_ref_np = np.array(U_X_REF, dtype=float)

strains_np = np.array(strains, dtype=float)
stresses_np = np.array(stresses, dtype=float)
displ_x_np = np.array(displacements, dtype=float)

print("\n--. stress-strain-displacement values --..")
print(f"{'Strain':<12} {'Stress (Pa)':<15} {'Ref Stress':<15} {'Displ (m)':<12} {'Plastic strain':<15}")
for e, s, sr, u, p in zip(strains, stresses, stresses_analytic_np, displacements, plastic_strains):
    print(f"{e:.4e}   {s:.4e}     {sr:.4e}     {u:.4e}     {p:.4e}")

# Plotting
plt.figure(figsize=(8, 6))
plt.plot(strain_ref_np, stresses_analytic_np, "k--", lw=1.5, label="Analytical (2D plane strain)")
plt.plot(strains_np, stresses_np, "r-o", lw=2, label="Numerical (Z3ST)")

# Mark yield point
E_eff = E / (1 - nu**2)
sigma_y_eff = sigma_y / np.sqrt(1 - nu + nu**2)
eps_y_eff = sigma_y_eff / E_eff
plt.plot(eps_y_eff, sigma_y_eff, 'ko', markersize=8, label="Analytic yield (2D)")

plt.xlabel(r"Strain $\epsilon_{xx}$ (/)")
plt.ylabel(r"Stress $\sigma_{xx}$ (/)")
plt.grid(True, linestyle=":", alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "stress_strain_curve.png"))
print("[INFO] stress_strain_curve.png saved\n")


# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---

# Strain
mask_zero = np.isclose(strain_ref_np, 0.0)
rel_error_E = np.zeros_like(strains_np)
rel_error_E[~mask_zero] = np.abs(strains_np[~mask_zero] - strain_ref_np[~mask_zero]) / np.abs(strain_ref_np[~mask_zero])
rel_error_E[mask_zero] = np.abs(strains_np[mask_zero] - strain_ref_np[mask_zero])

# Stress (compared against the analytical reference)
rel_error_S = np.zeros_like(stresses_np)
mask_zero_S = np.isclose(stresses_analytic_np, 0.0)
rel_error_S[~mask_zero_S] = np.abs(stresses_np[~mask_zero_S] - stresses_analytic_np[~mask_zero_S]) / np.abs(stresses_analytic_np[~mask_zero_S])
rel_error_S[mask_zero_S] = np.abs(stresses_np[mask_zero_S] - stresses_analytic_np[mask_zero_S])

# Displacement
rel_error_U = np.zeros_like(displ_x_np)
mask_zero_U = np.isclose(displ_x_ref_np, 0.0)
rel_error_U[~mask_zero_U] = np.abs(displ_x_np[~mask_zero_U] - displ_x_ref_np[~mask_zero_U]) / np.abs(displ_x_ref_np[~mask_zero_U])
rel_error_U[mask_zero_U] = np.abs(displ_x_np[mask_zero_U] - displ_x_ref_np[mask_zero_U])


errors = {
    "epsilon_xx": {
        "numerical": strains_np.tolist(),
        "reference": strain_ref_np.tolist(),
        "rel_error": rel_error_E.tolist(),
    },
    "sigma_xx": {
        "numerical": stresses_np.tolist(),
        "reference": stresses_analytic_np.tolist(), 
        "rel_error": rel_error_S.tolist(),
    },
    "u_xx": {
        "numerical": displ_x_np.tolist(),
        "reference": displ_x_ref_np.tolist(),
        "rel_error": rel_error_U.tolist(),
    },
}

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] non-regression completed.\n")
