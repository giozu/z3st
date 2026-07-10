"""
Calcul de l'évolution de la pression de contact Pk(t) pour un assemblage
fretté (arbre/moyeu) soumis au fluage, d'après Esposito, Bruno, Bertocco
(Int. J. Pressure Vessels and Piping 185 (2020) 104126).

Modèle : eq. (21)  ->  Pk(t) = Δu_el(t) * f
avec    Δu_el(t) = Δ * (1 + φ(1-n)*t / Δ^(1-n))^(1/(1-n))
et      φ = -(A/b) * f^n * (k1 + k2)
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import glob
import warnings

import yaml
import pyvista as pv

from z3st.utils.utils_load import generate_power_history

try:
    from z3st.utils.utils_extract_xdmf import *
except ModuleNotFoundError as exc:
    if exc.name != "h5py":
        raise
    extract_field_xdmf = None
    list_steps_xdmf = None


HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "output")
HISTORY_CSV = os.path.join(OUT, "history.csv")

# =====================================================================
# 1. VARIABLES A DEFINIR
# =====================================================================

R_PELLET = 0.0045
R_CLAD_I = 0.00453
R_CLAD_O = 0.005315
G0 = R_CLAD_I - R_PELLET
K_PEN = 5.0e13
T_TOTAL=1.5552e+8    # 1800 days
NB_STEPS=45
E_EL=99.3e9
NU=0.37
E_FUEL=200.0e9
NU_FUEL=0.345
ALPHA=7.5e-6
G=E_EL/(2*(1+NU))
CREEP_A0 = 2.82e-24 
CREEP_N = 3.0         # (-)
CREEP_Q = 1.2e5
SIGMA0=E_EL*5.0e-5/(R_CLAD_I-R_CLAD_O)
R_GAS=8.3145

# --- Géométrie du joint [mm] ---
a = 0.0      # rayon intérieur de l'arbre (a = 0 pour un arbre plein)
b = R_CLAD_I    # rayon commun (interface de contact)
c = R_CLAD_O     # rayon extérieur du moyeu

# --- Propriétés élastiques ---
E1, nu1 = E_FUEL, NU_FUEL   # arbre  [MPa, -]
E2, nu2 = E_EL, NU   # moyeu  [MPa, -]

# --- Propriétés de fluage (loi de Norton : eps_eq_dot = A * sigma_eq^n) ---
# On utilise ici les propriétés du moyeu (A2, n2). Si l'arbre est un
# arbre plein (a=0), k1=0 et les propriétés de l'arbre n'interviennent pas.
A1, n1 = 2e-14, 3.5     # arbre  [h^-1, -]
A2, n2 = 2e-14, 3.5     # moyeu  [h^-1, -]
n = CREEP_N                  # exposant utilisé dans la loi (n1 = n2 supposé)
A = CREEP_A0               # paramètre de Norton utilisé dans phi (relatif au moyeu)

# --- Tableau temps / interférence extrait du fichier history.csv ---
# t_array : temps [days]
# Delta   : interférence radiale [mm] (scalaire ou tableau)
def load_history_interference():
    """Charge le gap et l'interférence depuis history.csv."""
    if not os.path.exists(HISTORY_CSV):
        raise FileNotFoundError(f"{HISTORY_CSV} not found. Run the case first.")

    data = np.genfromtxt(HISTORY_CSV, delimiter=",", names=True, dtype=float)
    time_s = data["time_s"]
    time_h = time_s / 3600.0
    time_days = time_s / 86400.0
    gap_um = data["gap_um"]
    pressure = data["contact_pressure_MPa"]
    interference_um = np.maximum(0.0, -gap_um)  # gap négatif -> interférence
    return time_h[:NB_STEPS], time_days[:NB_STEPS], gap_um[:NB_STEPS], interference_um[:NB_STEPS], pressure[:NB_STEPS]


t_array, time_days, gap_um, interference_um, pressure = load_history_interference()
Delta = np.maximum(interference_um * 1e-6, 0.0)  # [mm] : 1 µm = 1e-3 mm


# =====================================================================
# 2. FONCTION DE CALCUL DE LA PRESSION DE CONTACT
# =====================================================================

def contact_pressure(t, Delta, a, b, c, E1, nu1, E2, nu2, A, n, T):
    """
    Calcule la pression de contact Pk(t) pour un joint fretté sous fluage.

    Paramètres
    ----------
    t : array_like
        Tableau des temps [h].
    Delta : float or array_like
        Interférence radiale [mm]. Peut être un scalaire ou un tableau.
    a, b, c : float
        Rayons intérieur, commun et extérieur [mm].
    E1, nu1 : float
        Module d'Young et coefficient de Poisson de l'arbre.
    E2, nu2 : float
        Module d'Young et coefficient de Poisson du moyeu.
    A, n : float
        Paramètres de la loi de Norton (moyeu).

    Retourne
    --------
    Pk : ndarray
        Pression de contact [MPa] pour chaque instant de t.
    f, k1, k2, phi : float
        Grandeurs intermédiaires du calcul.
    """
    t = np.asarray(t, dtype=float)
    Delta = np.asarray(Delta, dtype=float)
    # Ensure Young moduli are in MPa for the formula (E may be defined in Pa)
    E1 = float(E1) / 1e6
    E2 = float(E2) / 1e6
    Delta_arr = np.full_like(t, float(Delta[0]), dtype=float) if Delta.ndim == 0 else np.broadcast_to(Delta, t.shape).astype(float)
    Delta_safe = np.maximum(Delta_arr, 0)

    # --- Etape 1 : facteur élastique f (eq. 4-5) ---
    # Pk = Delta_u_el / [ b/E1*((a^2+b^2)/(b^2-a^2)+nu1) + b/E2*((c^2+b^2)/(c^2-b^2)-nu2) ]
    denom = (b / E1) * ((a**2 + b**2) / (b**2 - a**2) + nu1) \
            + (b / E2) * ((c**2 + b**2) / (c**2 - b**2) - nu2)
    f = 1.0 / denom  # [MPa/mm]

    # --- Etape 2 : k1, contribution de fluage de l'arbre (eq. 16) ---
    # k1 = 0 si arbre plein (a = 0) -> compression quasi-hydrostatique
    if a == 0:
        k1 = 0.0
    else:
        k1 = (np.sqrt(3) / 2) * (
            (np.sqrt(3) / n) * (a**(2/n) * b**(2/n)) / (b**(2/n) - a**(2/n))
        ) ** n

    # --- Etape 3 : k2, contribution de fluage du moyeu (eq. 20) ---
    k2=(np.sqrt(3) / 2) * (
            (np.sqrt(3) / n) * (c**(2/n) * b**(2/n)) / (c**(2/n) - b**(2/n))
        ) ** n

    # --- Etape 4 : phi, taux de relaxation (eq. 21-22) ---
    phi = -(CREEP_A0 * np.exp(-CREEP_Q / (R_GAS * T)) / b) * f**n * (k2 - k1)

    # --- Etape 5 : evolution de l'interference elastique Delta_u_el(t) (eq. 21) ---
    if n == 1:
        # cas particulier n=1 : pas de relaxation (cf. discussion de l'article)
        Delta_u_el = Delta_arr.copy()
    else:
        # Eviter la division par zéro lorsque Delta_safe == 0 (pas d'interférence)
        Delta_u_el = np.zeros_like(Delta_arr, dtype=float)
        positive = Delta_safe > 0
        if np.any(positive):
            # calcul uniquement pour les instants où l'interférence est strictement positive
            numer = 1 + phi[positive] * (1 - n) * t[positive] / (Delta_safe[positive] ** (1 - n))
            Delta_u_el[positive] = Delta_arr[positive] * (numer) ** (1 / (1 - n))

    # --- Etape 6 : pression de contact Pk(t) = Delta_u_el(t) * f (eq. 5/21) ---
    Pk = Delta_u_el * f
    return Pk, f, k1, k2, phi

def extract_clad_point_vonmises(xdmf_file,
                                r_target=(R_CLAD_I+R_CLAD_O)/2,
                                z_target=0.005):
    """
    Extract von Mises stress and temperature at one material point
    located in the middle of the cladding thickness.
    """

    if extract_field_xdmf is None:
        data = np.genfromtxt(HISTORY_CSV, delimiter=",", names=True, dtype=float)
        temp = np.asarray(data["T_max_K"], dtype=float)
        return np.zeros_like(temp), temp

    vm = []
    Temp = []

    try:
        steps = list_steps_xdmf(xdmf_file, "Stress")
    except Exception:
        steps = None

    n_steps = len(steps) if steps is not None else NB_STEPS + 1
    if steps is not None and n_steps == 0:
        raise ValueError(f"No saved steps found for 'Stress' in {xdmf_file}")

    try:
        for k in range(n_steps):

            r, z, _, S = extract_field_xdmf(
                xdmf_file,
                "Stress",
                step_index=k
            )

            xT, yT, zT, T = extract_field_xdmf(
                xdmf_file,
                "Temperature",
                step_index=k
            )

            S = np.asarray(S).reshape(-1,3,3)

            #
            # recherche du point le plus proche
            #
            d2 = (r-r_target)**2 + (z-z_target)**2

            i = np.argmin(d2)

            sigma = S[i]

            #
            # déviateur
            #
            p = np.trace(sigma)/3.0

            s = sigma-p*np.eye(3)


            sigma_eq = np.sqrt(1.5*np.sum(s*s))

            vm.append(sigma_eq)

            Temp.append(T[i])
    except (OSError, ValueError) as exc:
        # HDF5 file may be corrupted /unreadable (bad object header)
        # Fallback: read temperatures from history.csv and return zero vm
        data = np.genfromtxt(HISTORY_CSV, delimiter=",", names=True, dtype=float)
        temp = np.asarray(data["T_max_K"], dtype=float)
        temp = temp[:n_steps]
        vm = np.zeros_like(temp)
        return vm, temp

    return np.array(vm), np.array(Temp)



# =====================================================================
# 3. CALCUL ET AFFICHAGE
# =====================================================================

XDMF_FILE = os.path.join(OUT, "fields.xdmf")

vm,Temp=extract_clad_point_vonmises(XDMF_FILE, R_CLAD_I, 0.005)
""" Temp=np.array([500 for i in range(NB_STEPS)]) """

Pk, f, k1, k2, phi = contact_pressure(
    t_array, Delta, a, b, c, E1, nu1, E2, nu2, A, n, Temp
)



print(f"Facteur élastique f  = {f:.5f} MPa/mm")
print(f"k1                   = {k1:.5e}")
print(f"k2                   = {k2:.5e}")
print()
print("Tableau extrait de history.csv :")
print(f"{'t [days]':>10} | {'gap [um]':>10} | {'interference [um]':>18}")
print("-" * 45)
step = max(1, len(t_array) // 20)
for ti, gap_i, inter_i in zip(time_days[::step], gap_um[::step], interference_um[::step]):
    print(f"{ti:10.1f} | {gap_i:10.3f} | {inter_i:18.3f}")
print()
print(f"{'t [days]':>10} | {'Pk [MPa]':>10}")
print("-" * 25)
for ti, Pi in zip(time_days[::5], Pk[::5]):
    print(f"{ti:10.1f} | {Pi:10.3f}")

# --- Graphique Pk = f(t) ---
plt.figure(figsize=(7, 5))
plt.plot(time_days, Pk, "r-s", markersize=4, label=f"theoretical pressure for n = {n}")
plt.plot(time_days, pressure, "o", label="pressure from history.csv")
plt.xlabel("t [days]")
plt.ylabel("Pk [MPa]")
plt.title("Evolution de la pression de contact sous fluage pour k_pen=3,5")
plt.grid(True)
plt.legend()
plt.tight_layout()
output_path = os.path.join(OUT, "contact_pressure_evolution.png")
plt.savefig(output_path, dpi=150)
plt.close()
print(f"Figure written to {output_path}")
