#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Verification of the penalty contact pressure against the analytical Lame
interference-fit pressure. The inner cylinder is heated uniformly, so its free thermal
expansion is exactly u(b) = alpha_f (T - T_ref) b and the outer cylinder does not expand. 
The interference

    delta = alpha_f (T_pellet - T_ref) * b - g0

then gives the exact Lame shrink-fit pressure for a solid cylinder in a tube:

    p = delta / { b [ (1/E_c)((c^2+b^2)/(c^2-b^2) + nu_c) + (1/E_f)(1 - nu_f) ] }

(plane-stress form, consistent with the axially-free pellet)
"""

import os
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

T_TOTAL=2.16e8
E_CLAD=99.3e9
NU_CLAD=0.37
E_FUEL=200.0e9
NU_FUEL=0.345
ALPHA=7.5e-6
G=E_CLAD/(2*(1+NU_CLAD))
CREEP_A0 = 2.82e-24 
CREEP_N = 3.0         # (-)
CREEP_Q = 0
R_GAS=8.3145


files = sorted(glob.glob(os.path.join(OUT, "fields_*.vtu")))

geo = yaml.safe_load(open(os.path.join(CASE, "geometry.yaml")))
b = float(geo["outer_radius_1"])     # pellet outer / interface radius
bci = float(geo["inner_radius_2"])
c = float(geo["outer_radius_2"])

inp = yaml.safe_load(open(os.path.join(CASE, "input.yaml")))
cc = inp["models"]["contact"]
g0 = float(cc["initial_gap"])
k_pen = float(cc["penalty_stiffness"])
fuel = yaml.safe_load(open(os.path.join(CASE, inp["materials"]["cyl_1"])))
Ef, nuf, af, Trf = (float(fuel[k]) for k in ("E", "nu", "alpha", "T_ref"))
clad = yaml.safe_load(open(os.path.join(CASE, inp["materials"]["cyl_2"])))
Ec, nuc = float(clad["E"]), float(clad["nu"])

# Lame interference-fit compliance (plane stress): delta = p * b * comp
comp = (1.0 / Ec) * ((c**2 + bci**2) / (c**2 - bci**2) + nuc) + (1.0 / Ef) * (1.0 - nuf)

surf = lambda v, r, r0: v[np.abs(r - r0) < 2e-5].mean()
amean = lambda v, r, lo, hi: (lambda m: np.sum(v[m] * r[m]) / np.sum(r[m]))((r >= lo) & (r <= hi))

T_pellet, p_z3st, p_lame, interference, Delta = [], [], [], [], []
for f in files:
    m = pv.read(f)
    r = m.points[:, 0]
    ur = m.point_data["Displacement"][:, 0]
    T = m.point_data["Temperature"]

    Tp = amean(T, r, 0.0, b)                                  # uniform pellet temperature
    T_pellet.append(Tp)

    gap = g0 + surf(ur, r, bci) - surf(ur, r, b)
    interference.append(max(0,-gap))
    p_z3st.append(k_pen * max(0.0, -gap) / 1e6)              # MPa
    delta = af * (Tp - Trf) * b - g0                          # exact interference
    Delta.append(delta)
    p_lame.append((delta / (b * comp) / 1e6) if delta > 0 else 0.0)

T_pellet, p_z3st, p_lame = map(np.array, (T_pellet, p_z3st, p_lame))

plt.figure(figsize=(7, 5))
time_steps = np.linspace(0,T_TOTAL,13)
plt.plot(time_steps, p_lame, "k--", lw=1.5, label="Analytical Lame interference (exact)")
plt.plot(time_steps, p_z3st, "r-o", lw=2, label="Z3ST penalty contact")
plt.xlabel("pellet temperature (K)")
plt.ylabel("contact pressure (MPa)")
plt.title("Contact pressure verification: Z3ST vs Lame (uniform ΔT)")
plt.grid(True, ls=":", alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT, "contact_pressure_verification.png"), dpi=150)
print("[INFO] contact_pressure_verification.png saved")

mask = p_lame > 1.0       ### selects the value of the pressure not null to calculate the error 
""" for i in range(5):
    print("\n")

print(mask) """
if mask.any():
    rel = np.abs(p_z3st[mask] - p_lame[mask]) / p_lame[mask]
    print(k_pen, " penalty stiffness factor used")
    print(f"[INFO] closed-gap steps: {mask.sum()}, mean rel. error vs Lame = {rel.mean() * 100:.1f}%")
print("[INFO] non-regression completed.\n")


def plot_stress_profiles():
    """Radial profile of sigma_rr and sigma_theta at the last step,
    Z3ST vs the analytical Lame interference-fit solution.

      * pellet (0 <= r <= b), solid cylinder under external pressure p, axially
        free  ->  sigma_rr = sigma_theta = -p   (uniform)
      * clad (bci <= r <= c), tube with internal pressure p, free outer ->
            sigma_rr(r)    = p bci^2/(c^2-bci^2) (1 - c^2/r^2)
            sigma_theta(r) = p bci^2/(c^2-bci^2) (1 + c^2/r^2)

    sigma_rr is continuous, -p across the gap interface, and sigma_theta is
    compressive in the pellet but tensile in the clad.
    """
    if not files or not mask.any():
        print("[INFO] no closed-gap step: skipping stress profile")
        return

    i = int(np.argmax(p_lame))
    p = p_lame[i]                                   # MPa, exact interference pressure
    m = pv.read(files[i])

    if "Stress (cells)" in m.cell_data:
        coords = m.cell_centers().points
        s = np.asarray(m.cell_data["Stress (cells)"]).reshape(-1, 9)
    elif "Stress (points)" in m.point_data:
        coords = m.points
        s = np.asarray(m.point_data["Stress (points)"]).reshape(-1, 9)
    else:
        print("[INFO] no Stress field in VTU: skipping stress profile")
        return

    Lz = float(geo["Lz"])
    rr_c, zz_c = coords[:, 0], coords[:, 1]
    band = np.abs(zz_c - 0.5 * Lz) < 0.25 * Lz      # mid-height slice
    rb = rr_c[band]
    srr = s[band, 0] / 1e6                           # tensor order (r, theta, z)
    stt = s[band, 4] / 1e6
    o = np.argsort(rb)
    rb, srr, stt = rb[o], srr[o], stt[o]
    pellet = rb <= b + 1e-9
    cladm = rb >= bci - 1e-9

    # analytical curves at the same pressure p
    rp = np.linspace(0.0, b, 50)
    rcl = np.linspace(bci, c, 80)
    kk = p * bci**2 / (c**2 - bci**2)
    srr_clad = kk * (1.0 - c**2 / rcl**2)
    stt_clad = kk * (1.0 + c**2 / rcl**2)

    fig, ax = plt.subplots(figsize=(7.5, 5))
    # analytic (lines)
    ax.plot(rp * 1e3, np.full_like(rp, -p), color="C0", ls="--", lw=1.5,
            label=r"$\sigma_{rr}$ analytic")
    ax.plot(rcl * 1e3, srr_clad, color="C0", ls="--", lw=1.5)
    ax.plot(rp * 1e3, np.full_like(rp, -p), color="C3", ls=":", lw=1.8,
            label=r"$\sigma_{\theta\theta}$ analytic")
    ax.plot(rcl * 1e3, stt_clad, color="C3", ls=":", lw=1.8)
    # Z3ST (markers)
    ax.plot(rb[pellet] * 1e3, srr[pellet], "o", color="C0", ms=4,
            label=r"$\sigma_{rr}$ Z3ST")
    ax.plot(rb[cladm] * 1e3, srr[cladm], "o", color="C0", ms=4)
    ax.plot(rb[pellet] * 1e3, stt[pellet], "s", color="C3", ms=4,
            label=r"$\sigma_{\theta\theta}$ Z3ST")
    ax.plot(rb[cladm] * 1e3, stt[cladm], "s", color="C3", ms=4)

    ax.axvspan(b * 1e3, bci * 1e3, color="0.85", alpha=0.7)
    ax.axhline(0, color="grey", lw=0.8, ls="-")
    ax.text((b + bci) / 2 * 1e3, ax.get_ylim()[1] * 0.9, "gap",
            ha="center", fontsize=8, color="0.4")
    ax.set_xlabel("radius r (mm)")
    ax.set_ylabel("stress (MPa)")
    ax.set_title(f"Radial / hoop stress vs Lame (p = {p:.1f} MPa, mid-height)")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(True, ls=":", alpha=0.5)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "stress_profile_verification.png"), dpi=150)
    plt.close(fig)
    print("[INFO] stress_profile_verification.png saved")

    # interface continuity diagnostic
    if pellet.any() and cladm.any():
        srr_fuel_surf = srr[pellet][-1]
        srr_clad_inner = srr[cladm][0]
        print(f"[INFO] interface sigma_rr: pellet={srr_fuel_surf:8.2f} MPa, "
              f"clad={srr_clad_inner:8.2f} MPa, analytic -p={-p:8.2f} MPa")
        print(f"[INFO] pellet sigma_rr range: [{srr[pellet].min():.2f}, "
              f"{srr[pellet].max():.2f}] MPa (should be ~ -p, uniform, never tensile)")


plot_stress_profiles()


def plot_interference_pressure_time():
    """Plot interference and p_z3st as functions of time (file index)."""
    time_steps = np.arange(len(interference))
    print(interference)

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(7, 9), sharex=True)
    
    # Interference plot
    ax1.plot(time_steps, interference, "b-o", lw=2, ms=4, label="Interference")
    ax1.set_ylabel("interference (m)")
    ax1.set_title("Interference and contact pressure evolution")
    ax1.grid(True, ls=":", alpha=0.6)
    ax1.legend()
    
    # Contact pressure plot
    ax2.plot(time_steps, p_z3st, "r-s", lw=2, ms=4, label="p_z3st")
    ax2.set_xlabel("time step")
    ax2.set_ylabel("contact pressure (MPa)")
    ax2.grid(True, ls=":", alpha=0.6)
    ax2.legend()

    ax3.plot(interference, p_z3st, "b-o", lw=2, ms=4, label="p_z3st vs Interference")
    ax3.set_ylabel("pressure contact")
    ax3.set_title("Interference and contact pressure evolution")
    ax3.grid(True, ls=":", alpha=0.6)
    ax3.legend()
    
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "interference_pressure_time.png"), dpi=150)
    plt.close(fig)
    print("[INFO] interference_pressure_time.png saved")


plot_interference_pressure_time()

Kpen=[1e12,5e12,1e13,5e13,1.1e14,1.2e14,1.3e14,1.4e14,1.5e14,2e14,5e14]

error_pressure=[7.45e-01,3.73e-01,2.34e-01,6.85e-02,4.35e-02,4.87e-02,3.52e-02,3.36e-02,0.19,1.74e+01,5.66e+01]
plt.figure(figsize=(8, 5))
plt.semilogx(Kpen, error_pressure, "o-", color="blue", lw=2, markersize=8, label="Erreur pression")
plt.axhline(0, color="black", linestyle="--", lw=1, label="Solution Lamé (erreur = 0)")
plt.xlabel("penalty stiffness $K_{pen}$ (Pa)")
plt.ylabel("Relative error on the contact pressure")
plt.title("Error vs penalty stiffness (gap = 65 µm)")
plt.grid(True, which="both", linestyle="--", alpha=0.5)
plt.legend()
plt.tight_layout()
plt.savefig("penalty_stiffness_error.png", dpi=150)
plt.show()

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
    E1 = float(E1) 
    E2 = float(E2)
    Delta_arr = np.full_like(t, float(Delta[0]), dtype=float) if Delta.ndim == 0 else np.broadcast_to(Delta, t.shape).astype(float)
    Delta_safe = np.maximum(Delta_arr, 0)

    # --- Etape 1 : facteur élastique f (eq. 4-5) ---
    # Pk = Delta_u_el / [ b/E1*((a^2+b^2)/(b^2-a^2)+nu1) + b/E2*((c^2+b^2)/(c^2-b^2)-nu2) ]
    denom = (b / E1) * ((a**2 + b**2) / (b**2 - a**2) + nu1) \
            + (b / E2) * ((c**2 + b**2) / (c**2 - b**2) - nu2)
    f = 1.0 / denom  # [MPa/mm]
    print(a,b,c)

    # --- Etape 2 : k1, contribution de fluage de l'arbre (eq. 16) ---
    # k1 = 0 si arbre plein (a = 0) -> compression quasi-hydrostatique
    if a == 0:
        k1 = 0.0
    else:
        k1 = (np.sqrt(3) / 2) * (
            (np.sqrt(3) / n) * (a**(2/n) * b**(2/n)) / (b**(2/n) - a**(2/n))
        ) ** n

    # --- Etape 3 : k2, contribution de fluage du moyeu (eq. 20) ---
    k2= (1- n*b**(2/n)/4 * ( (c**(2/n) - b**(2/n))/ c**(2/n) * b**(2/n))) * ( 
            (2 / n) * (c**(2/n) * b**(2/n)) / (c**(2/n) - b**(2/n))
        ) ** n

    # --- Etape 4 : phi, taux de relaxation (eq. 21-22) ---

    phi = -(CREEP_A0 * np.exp(-CREEP_Q / (R_GAS * T)) / b) * f**n * (k2 + k1)
    print(f,k1,k2)
    contact_mask = Delta_arr > 1e-9          # seuil physique, pas 1e-12
    Pk = np.zeros_like(t)

    if np.any(contact_mask):
        idx0 = np.argmax(contact_mask)        # premier indice de contact
        Delta0 = Delta_arr[idx0]              # interférence "initiale" du fluage
        t0 = t[idx0]                          # instant où le fluage démarre
        t_rel = t[idx0:] - t0   

    # --- Etape 5 : evolution de l'interference elastique Delta_u_el(t) (eq. 21) ---
    if n == 1:
        # cas particulier n=1 : pas de relaxation (cf. discussion de l'article)
        Delta_u_el = Delta_arr.copy()
    else:
        # Eviter la division par zéro lorsque Delta_safe == 0 (pas d'interférence)
        Delta_u_el = np.zeros_like(Delta_arr, dtype=float)

        numer = 1 + phi[idx0:] * (1 - n) * t_rel / (Delta_safe[idx0:] ** (1 - n))
        Delta_u_el[idx0:] = Delta_arr[idx0:] * (numer) ** (1 / (1 - n))

    # --- Etape 6 : pression de contact Pk(t) = Delta_u_el(t) * f (eq. 5/21) ---
    Pk = Delta_u_el * f
    return Pk/1e6, f, k1, k2, phi


Temp=np.array([600 for i in range(len(time_steps))])
P_c, f, k1, k2, phi = contact_pressure(time_steps,Delta,0,bci,c,E_FUEL,NU_FUEL,E_CLAD, NU_CLAD, CREEP_A0, CREEP_N, Temp)


print(Delta)
print(P_c)
print(time_steps)
plt.figure(figsize=(7, 5))
plt.plot(time_steps, P_c, "k--", lw=1.5, label="Analytical Lame interference (exact)")
plt.plot(time_steps, p_z3st, "r-o", lw=2, label="Z3ST penalty contact")
plt.xlabel("Time ( days )")
plt.ylabel("contact pressure (MPa)")
plt.title("Contact pressure verification: Z3ST vs Eposito's solution (uniform ΔT)")
plt.grid(True, ls=":", alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT, "pressure_creep_analysis.png"), dpi=150)
print("[INFO] contact_pressure_verification.png saved")




# --. numerical results --..
errors = {
    "contact_pressure": {
        "numerical": p_z3st[mask].max(),
        "reference": p_lame[mask].max(),
        "abs_error": float(rel.max()),
        "rel_error": float(rel.max()),
    },
}

pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE)
regression_check(errors, CASE)

print("\n[INFO] non-regression completed.\n")
