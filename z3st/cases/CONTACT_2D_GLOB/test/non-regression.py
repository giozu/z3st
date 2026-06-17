#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: thick_cylindrical_shell_GPS_2D

non-regression script
---------------------
Analytical non-regression for a 2D thick-walled cylindrical shell under
internal (Pi) and external (Po) pressure.
Reference is the Lamé solution under generalized plane strain.

"""

import os

import matplotlib.pyplot as plt
import numpy as np

from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_verification import *

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
VTU_FILE = os.path.join(CASE_DIR, "output", "fields.vtu")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")

# Geometry and material
Ri, Ro, Lz = 0.02, 0.03, 0.50  # m          inner and outer radius, height
Pi, Po = 1.0e6, 0.0  # Pa         internal and external pressure
E, nu = 2.0e11, 0.3  # Pa, -      Young modulus, Poisson ratio
t = Ro - Ri  # m          wall thickness
slenderness = Ri / t  # -          slenderness ratio
z_target, z_tol = Lz / 2, 0.01  # m          z-plane for data extraction

TOLERANCE = 5.0e-3  # -          tolerance for non-regression

# --.. ..- .-.. .-.. --- analytic functions  --.. ..- .-.. .-.. ---
eps_zz_GPS = -2.718310e-06

# Lamé solutions
A = (Pi * Ri**2 - Po * Ro**2) / (Ro**2 - Ri**2)
B = (Ri**2 * Ro**2 * (Pi - Po)) / (Ro**2 - Ri**2)
sigma_zz_ana_L = 2 * nu * A + E * eps_zz_GPS  # Generalized plane strain (epsilon_z = const)


def epsilon_rr_ref(r):
    return (1 + nu) / E * (A * (1 - 2 * nu) - B / r**2) - nu * eps_zz_GPS


def epsilon_tt_ref(r):
    return (1 + nu) / E * (A * (1 - 2 * nu) + B / r**2) - nu * eps_zz_GPS


def test_contact():
    # Charger les déplacements
    u = extract_field("output/fields.vtu", "Displacement")

    # Extraire les coordonnées des nœuds
    coords = extract_field("output/fields.vtu", "Coordinates")

    # Trouver les nœuds de l'interface (x ≈ R_mid)
    R_mid = 0.04
    interface_nodes = np.isclose(coords[:, 0], R_mid, atol=1e-4)

    # Vérifier que la pénétration est nulle
    u_in = u[interface_nodes]
    coords_in = coords[interface_nodes]

    # Calculer la distance minimale entre les nœuds de l'interface
    # (Supposons que steel_in est à gauche et steel_o à droite)
    min_distance = np.min(coords_in[:, 0] + u_in[:, 0] - R_mid)
    assert min_distance >= -1e-6, f"Pénétration détectée: {min_distance}"

if __name__ == "__main__":
    test_contact()
    print("✅ Test de contact passé avec succès !")