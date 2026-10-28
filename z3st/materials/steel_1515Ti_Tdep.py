# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: temperature-dependent properties of 15-15Ti steel (20 % CW)
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
15-15Ti stabilised austenitic steel, 20 % cold worked: property correlations
for the card ``materials/1515Ti_Tdep.yaml``.

Two correlation sets are given for most properties, [Sobolev] and
[DiMarcello]. T is in K everywhere in the signatures; the °C correlations
convert internally (t = T - 273.15).

Hooks the card can point to (resolved by ``spine.load_materials``):

- ``k``: :func:`k_Sobolev`, :func:`k_DiMarcello` (W/(m·K)); lagged on the
  staggered iterate T_iter, as every string k card.
- ``cp``: :func:`cp_Sobolev`, :func:`cp_DiMarcello` (J/(kg·K)), UFL(T_iter).
- ``E``: :func:`E_Sobolev`, :func:`E_DiMarcello` (Pa), UFL(T_iter).
- ``eigenstrain``: :func:`eigenstrain`, thermal expansion η(T) and,
  optionally, void swelling at a uniform fluence. Replaces ``alpha``: a card
  using it must not carry ``alpha``, or the thermal strain is counted twice.
- ``creep: nabarro_herring_structure`` with ``creep_A0``, ``creep_B``, ``creep_Q`` in SI from
  :func:`sinh_creep_SI` (values for 1.4970 and AIM1 in the card).

The other functions (ρ, YS, UTS, creep rates, rupture time, CDF) are for
post-processing: they accept floats or numpy arrays, the polynomial ones and
the creep rates also UFL. ρ(T) never enters the solver: the thermal problem is written on the
undeformed mesh, where the mass per unit volume is the reference ρ0 (card
``rho``); ρ(T) would count the thermal expansion twice in ρ·cp.

Card keys read by :func:`eigenstrain`
-------------------------------------
- ``T_ref``          (K)  reference temperature of η, 293.15 (t_ref = 20 °C)
- ``expansion_law``      ``Sobolev`` (default) | ``DiMarcello`` | ``none``
- ``swelling_law``       ``none`` (default) | ``generalized`` | ``aim1``
- ``fluence``            (1e22 n/cm²) uniform fast fluence, for the swelling
- ``axial_creep_law``    ``none`` (default) | ``1.4970`` | ``AIM1``
- ``axial_creep_sigma``  (MPa) imposed uniaxial stress
- ``axial_creep_time``   (h) time at that stress
- ``axial_creep_axis``   index of the stress axis, 0 | 1 | 2 (default 2, z in 3d)

Each contribution is switched off on its own (``none``), so thermal bowing,
swelling and creep can be run separately and then together.
Imposed creep vs ``creep: nabarro_herring_structure``: the imposed creep
takes σ from the card and never relaxes it; ``creep:
nabarro_herring_structure`` (models/creep_model.py) takes σ_eq from the FE
solution and relaxes it. Do not use both on one card.

Sources
-------
[Sobolev]    V. Sobolev, as collected in the course notes "1515Ti stabilized
             20 % cold worked" (2025).
[DiMarcello] V. Di Marcello et al. (TRANSURANUS 15-15Ti), same notes.
[Tobbe]      H. Többe, 1975, 1.4970 creep and rupture, same notes.
[Banerjee]   A. Banerjee, S. Raju, R. Divakar, E. Mohandas, High Temperature
             Heat Capacity of Alloy D9 Using Drop Calorimetry Based Enthalpy
             Increment Measurements, Int. J. Thermophys. 28 (2007) 97-108,
             doi:10.1007/s10765-006-0136-0.
[Thesis]     "Extension and application of a fuel pin performance code in
             the conceptual design of a Lead-cooled Fast Reactor", Master
             Thesis in Nuclear Engineering, Politecnico di Milano:
             Tables 1.5 (η), 1.6 (cp), 1.7 (λ), 1.8 (E), 1.9 (YS, UTS),
             eq. 1.9 (YS), 1.11-1.12 (thermal creep), 1.14 (irradiation
             creep); C = 20 of the rupture parameter.
[own]        Derivation or correction written for this module.

Corrections to [Thesis], marked [TBC] where they occur
------------------------------------------------------
- k [DiMarcello]: Table 1.7 gives (202700 - 81.67 T)/1000, the E law of
  Table 1.8 (139.6 W/(m·K) at 773 K). Replaced by 13.95 + 0.01163 t, the
  TRANSURANUS law already quoted in ``15_15Ti.yaml``.
"""

import numpy as np
import ufl

T_C0 = 273.15                 # (K) 0 °C
RHO0 = 7900.0                 # (kg/m³) density at t_ref = 20 °C
R_CAL = 1.986                 # (cal/(mol·K)) gas constant of the creep laws
R_GAS = 8.314462618           # (J/(mol·K)) gas constant, as in models/creep_model.py
FLUENCE_PER_DPA = 2.0e21      # (n/cm²) per dpa
FLUENCE_RUPTURE = 1.0e22      # (n/cm²) above it the rupture strain is 6 %
RUPTURE_STRAIN_IRR = 0.06     # (-)


def _is_ufl(x):
    return isinstance(x, ufl.core.expr.Expr)


def _exp(x):
    return ufl.exp(x) if _is_ufl(x) else np.exp(x)


def _sinh(x):
    return ufl.sinh(x) if _is_ufl(x) else np.sinh(x)


def _num(x):
    """UFL expressions pass through; numbers and lists become float arrays."""
    return x if _is_ufl(x) else np.asarray(x, dtype=float)


# --. thermal --..


def k_Sobolev(T):
    """λ(T) (W/(m·K)) [Sobolev]: 8.826 + 1.707e-2 T - 2.315e-6 T².

    Source: [Thesis] Table 1.7 (Sobolev 2011; Rudtsch et al. 2005).
    """
    return 8.826 + 1.707e-2 * T - 2.315e-6 * T**2


def k_DiMarcello(T):
    """λ(T) (W/(m·K)) [DiMarcello]: 13.95 + 0.01163 t, t in °C.

    [TBC] [Thesis] Table 1.7 repeats the E law here; see the module docstring.
    """
    return 13.95 + 0.01163 * (T - T_C0)


def cp_Sobolev(T):
    """cp(T) (J/(kg·K)) [Sobolev]: 393.3 + 0.3609 T - 2.431e-4 T² + 9.202e-8 T³.

    Source: [Thesis] Table 1.6 (Sobolev 2011; [Banerjee]).
    """
    return 393.3 + 0.3609 * T - 2.431e-4 * T**2 + 9.202e-8 * T**3


def cp_DiMarcello(T):
    """cp(T) (J/(kg·K)) [DiMarcello]: 431 + 0.177 T + 8.72e-5 T^-2, T in K.

    Source: [Thesis] Table 1.6 (Di Marcello 2013), as printed. The T^-2 term
    is ~1e-10 J/(kg·K) in any T range of interest.
    """
    return 431.0 + 0.177 * T + 8.72e-5 * T**-2


# --. thermal expansion and density --..


def eta_Sobolev(T):
    """η = ΔL/L0 (-) [Sobolev], T in K:
    (-4.019 + 1.105e-2 T + 1.191e-5 T² - 8.879e-9 T³ + 3.242e-12 T⁴)·1e-3.

    Source: [Thesis] Table 1.5 (Sobolev 2011; Nikolopoulos and Schulz 1979).
    """
    return (-4.019 + 1.105e-2 * T + 1.191e-5 * T**2
            - 8.879e-9 * T**3 + 3.242e-12 * T**4) * 1.0e-3


def eta_DiMarcello(T):
    """η = ΔL/L0 (-) [DiMarcello]: -3.101e-4 + 1.545e-5 t + 2.75e-9 t², t in °C."""
    t = T - T_C0
    return -3.101e-4 + 1.545e-5 * t + 2.75e-9 * t**2


_ETA = {"sobolev": eta_Sobolev, "dimarcello": eta_DiMarcello}


def _eta(law):
    return _ETA[str(law).lower()]


def rho(T, law="Sobolev"):
    """ρ(T) = ρ0 / (1 + η(T))³ (kg/m³), ρ0 = 7900 kg/m³ at t_ref = 20 °C.

    Mass conservation in a cube whose side scales by (1 + η), η from 20 °C.
    Post-processing only (see the module docstring): the card's ``rho`` is ρ0.
    """
    return RHO0 / (1.0 + _eta(law)(T)) ** 3


def thermal_strain(T, T_ref, law="Sobolev"):
    """Linear thermal strain from T_ref, (η(T) - η(T_ref)) / (1 + η(T_ref)) (-).

    Source: [own], the length change referred to the length at T_ref.
    """
    eta = _eta(law)
    eta_ref = float(eta(T_ref))
    return (eta(T) - eta_ref) / (1.0 + eta_ref)


# --. elasticity and strength --..


def E_Sobolev(T):
    """E(T) (Pa) [Sobolev]: (213.8 - 5.545e-2 T - 2.039e-5 T²) GPa.

    Source: [Thesis] Table 1.8 (Sobolev 2011; Bergmann 1988).
    """
    return (213.8 - 5.545e-2 * T - 2.039e-5 * T**2) * 1.0e9


def E_DiMarcello(T):
    """E(T) (Pa) [DiMarcello]: (202700 - 81.67 t) / 1000 GPa, t in °C.

    Source: [Thesis] Table 1.8 (Di Marcello 2013).
    """
    return (202700.0 - 81.67 * (T - T_C0)) * 1.0e6


def yield_strength_Sobolev(T):
    """YS(T) (MPa) [Sobolev]:
    1560 - 5.86 T + 0.012 T² - 1.04e-5 T³ + 3.06e-9 T⁴.

    Source: [Thesis] Table 1.9 (Sobolev 2011; Bergmann 1988).
    """
    return 1560.0 - 5.86 * T + 0.012 * T**2 - 1.04e-5 * T**3 + 3.06e-9 * T**4


def uts_Sobolev(T):
    """UTS(T) (MPa) [Sobolev]:
    1410 - 4.35 T + 9.03e-3 T² - 8.14e-6 T³ + 2.42e-9 T⁴.

    Source: [Thesis] Table 1.9 (Sobolev 2011; Bergmann 1988).
    """
    return 1410.0 - 4.35 * T + 9.03e-3 * T**2 - 8.14e-6 * T**3 + 2.42e-9 * T**4


def yield_strength_DiMarcello(T):
    """YS(T) (MPa) [DiMarcello], t in °C:
    555.5 - 0.25 t (t < 600), 405.5 - 0.775 (t - 600) (600 ≤ t < 1000),
    345.5 - 0.25 t (t ≥ 1000).

    Source: [Thesis] eq. 1.9. Continuous at 600 °C and 1000 °C.
    """
    t = np.asarray(T, dtype=float) - T_C0
    return np.where(t < 600.0, 555.5 - 0.25 * t,
                    np.where(t < 1000.0, 405.5 - 0.775 * (t - 600.0),
                             345.5 - 0.25 * t))


def rupture_strain_irradiated(fluence):
    """Rupture strain (-) in pile: 6 % once the fast fluence exceeds 1e22 n/cm²,
    NaN below (the notes give no value there). ``fluence`` in n/cm².
    """
    return np.where(np.asarray(fluence) > FLUENCE_RUPTURE, RUPTURE_STRAIN_IRR, np.nan)


# --. void swelling --..


def fluence_from_dpa(dpa):
    """Fast fluence Φ (1e22 n/cm²) from the dose (dpa), 2e21 n/cm² = 1 dpa."""
    return dpa * FLUENCE_PER_DPA / 1.0e22


def swelling_generalized(T, fluence):
    """ΔV/V (-) "Generalized 15-15Ti" (older steel), Φ in 1e22 n/cm²:
    1.5e-3·exp(-2.5((t - 450)/100)²)·Φ^2.75 %.
    """
    t = T - T_C0
    return 1.5e-3 * _exp(-2.5 * ((t - 450.0) / 100.0) ** 2) * fluence**2.75 / 100.0


def swelling_aim1(T, fluence):
    """ΔV/V (-) "Specific AIM1" (Phénix), Φ in 1e22 n/cm²:
    1.3e-5·exp(-((t - 490)/100)²)·Φ^3.9 %.
    """
    t = T - T_C0
    return 1.3e-5 * _exp(-(((t - 490.0) / 100.0) ** 2)) * fluence**3.9 / 100.0


_SWELLING = {"generalized": swelling_generalized, "aim1": swelling_aim1}


def axial_creep_strain(T, material, dim=3):
    """Creep strain (-) under an imposed uniaxial stress σ along axis a,

        ε_cr = ε̇(T, σ)·t·(e_a⊗e_a - ½ (I - e_a⊗e_a)),

    i.e. ε̇·t along a and -½ ε̇·t across it (volume-preserving flow, the
    uniaxial case of von Mises flow, σ_eq = σ). ε̇ (1/s) from
    :func:`creep_rate_1_4970` or :func:`creep_rate_aim1` with
    σ = ``axial_creep_sigma`` (MPa); t = ``axial_creep_time`` (h) · 3600 s.

    Source: [own] closed form of ∫ ε̇ dt at constant T and σ.
    Validity: steady-state T (exact; with T varying in time it would need
    the T history) and σ fixed by the card, not updated by the FE stress.
    """
    law = str(material.get("axial_creep_law", "none")).lower()
    sigma = float(material["axial_creep_sigma"])
    t_s = float(material["axial_creep_time"]) * 3600.0
    a = int(material.get("axial_creep_axis", 2))
    eps_cr = _CREEP_RATE[law](T, sigma) * t_s
    e = ufl.unit_vector(a, dim)
    P = ufl.outer(e, e)
    return eps_cr * (P - 0.5 * (ufl.Identity(dim) - P))


def eigenstrain(T, material, model=None, dim=3):
    """Eigenstrain ε* (-) = ε_th(T) I + ((1 + ΔV/V)^(1/3) - 1) I + ε_cr.

    Card hook: ``eigenstrain: materials.steel_1515Ti_Tdep.eigenstrain``, in
    place of ``alpha``. Three contributions, each switched off by ``none``:
    - ε_th from :func:`thermal_strain` with ``expansion_law``;
    - ΔV/V from ``swelling_law`` at the uniform card ``fluence``
      (1e22 n/cm²);
    - ε_cr from :func:`axial_creep_strain` with ``axial_creep_law``.
    All are evaluated on the current T (path-independent, exact at
    steady-state T).
    Validity: small strains, isotropic expansion and swelling, no stress
    dependence of swelling. With the thermal model off (T is None) only the
    imposed creep remains, at the card ``T_initial``.
    """
    law = str(material.get("expansion_law", "Sobolev"))
    if "alpha" in material and law.lower() != "none":
        raise ValueError("steel_1515Ti_Tdep.eigenstrain includes η(T): remove 'alpha' from the card "
                         "or set expansion_law: none.")
    I = ufl.Identity(dim)
    eps = 0.0 * I
    if str(material.get("axial_creep_law", "none")).lower() != "none":
        T_cr = T if T is not None else float(material["T_initial"])
        eps = eps + axial_creep_strain(T_cr, material, dim)
    if T is None:
        return eps

    if law.lower() != "none":
        eps = eps + thermal_strain(T, float(material["T_ref"]), law) * I

    sw = str(material.get("swelling_law", "none")).lower()
    if sw != "none":
        g = _SWELLING[sw](T, float(material["fluence"]))
        eps = eps + ((1.0 + g) ** (1.0 / 3.0) - 1.0) * I
    return eps


# --. thermal creep and rupture --..


def sinh_creep_SI(A_pct_h, Q_cal, B_cal):
    """Card values (creep_A0, creep_Q, creep_B) of ``creep: nabarro_herring_structure`` from a
    thermal-creep law in the units of [Thesis] Tables 1.11-1.12:

        ε̇ (%/h) = A·exp(-Q/(R_c·T))·sinh(B·σ/(R_c·T)),
        σ in MPa, T in K, Q in cal/mol, B in cal/(mol·MPa), R_c = 1.986 cal/(mol·K).

    z3st (models/creep_model.py) evaluates, σ in Pa, R = R_GAS in J/(mol·K):

        ε̇ (1/s) = creep_A0·exp(-creep_Q/(R·T))·sinh(creep_B·σ/(R·T)).

    Conversions:
    - creep_A0 (1/s)         = A / 100 / 3600            (%/h → 1/s)
    - creep_Q  (J/mol)       = Q · R_GAS/R_c             (keeps Q/R in K)
    - creep_B  (J/(mol·Pa))  = B · R_GAS/R_c · 1e-6      (keeps B/R; MPa → Pa)

    [TBC] rule 4 — cut to "Q and B are scaled by R_GAS/R_c = 4.1865 J/cal."?
    loses: the warning not to use 4.184 J/cal and its ≈ 3 % effect on ε̇.
    Q and B are scaled by R_GAS/R_c = 4.1865 J/cal rather than 4.184, so the
    exponents Q/(R·T) and B·σ/(R·T) are those of the source: the 0.06 %
    difference, times Q/(R·T) ≈ 50, would change ε̇ by ≈ 3 %.
    """
    scale = R_GAS / R_CAL
    return A_pct_h / 100.0 / 3600.0, Q_cal * scale, B_cal * scale * 1.0e-6



def creep_rate_1_4970(T, sigma):
    """Thermal creep rate (1/s) of 1.4970 [Tobbe], σ in MPa, T in K:
    7.49e15·exp(-91000/(R T))·sinh(57.46 σ/(R T)) %/h, R = 1.986 cal/(mol·K),
    Q = 91000 cal/mol, B = 57.46 cal/(mol·MPa); returned as %/h / 100 / 3600.

    Source: [Thesis] Table 1.11 (Többe 1975, TRANSURANUS v1m3j12).
    """
    RT = R_CAL * _num(T)
    pct_h = 7.49e15 * _exp(-91000.0 / RT) * _sinh(57.46 * _num(sigma) / RT)
    return pct_h / 100.0 / 3600.0


def creep_rate_aim1(T, sigma):
    """Thermal creep rate (1/s) of AIM1 (fit), σ in MPa, T in K:
    2.30e14·exp(-84600/(R T))·sinh(39.72 σ/(R T)) %/h, R = 1.986 cal/(mol·K),
    Q = 84600 cal/mol, B = 39.72 cal/(mol·MPa); returned as %/h / 100 / 3600.

    Source: [Thesis] Table 1.12 (fit of Gavoille et al. 2013 data).
    """
    RT = R_CAL * _num(T)
    pct_h = 2.30e14 * _exp(-84600.0 / RT) * _sinh(39.72 * _num(sigma) / RT)
    return pct_h / 100.0 / 3600.0


# "1.497": an unquoted 1.4970 in YAML is read as the float 1.497
_CREEP_RATE = {"1.4970": creep_rate_1_4970, "1.497": creep_rate_1_4970,
               "aim1": creep_rate_aim1}


def rupture_parameter_1_4970(sigma):
    """Rupture parameter P (K) of 1.4970 from σ (MPa), Nabarro-Herring form:
    0.01892153 σ² - 24.19134 σ + 21761.81 (σ < 395),
    -6.017016 σ + 17535.01 (σ ≥ 395).
    """
    s = np.asarray(sigma, dtype=float)
    return np.where(s < 395.0,
                    0.01892153 * s**2 - 24.19134 * s + 21761.81,
                    -6.017016 * s + 17535.01)


def rupture_time_1_4970(T, sigma, C=20.0):
    """Rupture time t_R (h) of 1.4970 from P = T (C + log10(t_R / T)).

    Source: C = 20 from [Thesis].
    """
    T = np.asarray(T, dtype=float)
    return T * 10.0 ** (rupture_parameter_1_4970(sigma) / T - C)


def rupture_time_aim1(T, sigma):
    """Rupture time t_R (h) of AIM1 from the Larson-Miller parameter,
    T (17.125 + log10 t_R) = (2060 - σ) / 0.095, σ in MPa, T in K.
    """
    T = np.asarray(T, dtype=float)
    return 10.0 ** ((2060.0 - np.asarray(sigma)) / (0.095 * T) - 17.125)


def cdf_time(dt_h, t_rupture_h):
    """Time-fraction cumulative damage CDF = Σ dt_i / t_R,i (-), dt in h."""
    return np.cumsum(np.asarray(dt_h) / np.asarray(t_rupture_h))


def cdf_strain(d_eps, eps_rupture):
    """Strain-fraction cumulative damage CDF = Σ Δε_i / ε_R,i (-)."""
    return np.cumsum(np.asarray(d_eps) / np.asarray(eps_rupture))


# --. irradiation creep --..


def creep_irr_B(E_avg):
    """Irradiation-creep coefficient B (Pa⁻¹ per n/m²) for the card key
    ``creep_irr_B`` (with ``fast_flux`` in n/(m²·s)), from
    ε̇ = 3.2e-24·Ē·φ·σ, Ē in MeV, φ in n/(cm²·s), σ in MPa, ε̇ in %/h.

    B = 3.2e-24·Ē·1e-4·1e-6 / 100 / 3600 = 8.89e-40·Ē.

    Source: [Thesis] Table 1.14 (Többe 1975, TRANSURANUS v1m3j12).
    With 2e21 n/cm² = 1 dpa and Ē = 0.5 MeV this is 0.9e-8 1/(MPa·dpa),
    100 times below the 1.1e-6 1/(MPa·dpa) of Ye et al. 2020 used in
    ``15_15Ti.yaml``.
    """
    return 3.2e-24 * E_avg * 1.0e-4 * 1.0e-6 / 100.0 / 3600.0
