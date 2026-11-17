# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.3.2 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Void swelling and irradiation creep of fast-reactor structural steels.

Both contributions are driven by the displacement-damage dose ``model.dose``
(dpa) and its rate ``model.dose_rate`` (dpa/s). ``spine.initialize_fields``
builds them for every material whose card carries ``dpa_rate``;
``spine.update_state`` accumulates the dose once per step, like the burnup of a
fissile material.

- Void swelling: an eigenstrain, ``eigenstrain:
  materials.steel_irradiation.void_swelling``. Depends on T and D.
- Irradiation creep: a creep compliance for CreepModel (with ``creep: norton``),
  chosen in the card by one of two hooks:

  * ``creep_irr: materials.steel_irradiation.irradiation_creep_uncoupled``
    ε̇ = B0·Ḋ·σ_eq, independent of swelling. Used in the first phase, where
    swelling and creep are kept separate;
  * ``creep_irr: materials.steel_irradiation.irradiation_creep``
    ε̇ = (B0·Ḋ + D_sc·ġ)·σ_eq, with the swelling-driven term.

Either contribution can be switched on alone.

Card keys (SI, T in K, dose in dpa)
-----------------------------------
- ``dpa_rate``         (dpa/s) dose rate (loading, set by the case)
- ``swelling_g0_max``  (1/dpa) peak steady-state volumetric swelling rate
- ``swelling_T_peak``  (K)     temperature of the peak rate
- ``swelling_T_width`` (K)     Gaussian width of the bell in T
- ``swelling_D_inc``   (dpa)   incubation dose
- ``swelling_D_tr``    (dpa)   transition dose; 0 → bilinear onset
- ``creep_irr_B0``     (1/(Pa·dpa)) irradiation-creep compliance
- ``creep_irr_Dsc``    (1/Pa)  creep-swelling coupling coefficient
                               (``irradiation_creep`` only)

Sources 
------------------------------------------------------
[Acosta2019] F. Acosta et al., Nucl. Eng. Des. 348 (2019) 90-106,
             Eq. 9-15 (DOMAJEUR2 cladding laws, CEA).
[Acosta2021] F. Acosta et al., Nucl. Eng. Des. 374 (2021) 111048,
             Eq. 1-3.
[Ye2020]     X. Ye et al., Nucl. Eng. Des. 364 (2020) 110676, Eq. 2-3
             (316Ti CW cladding/wire/duct; test subassembly at 560 °C,
             peak 75 dpa, 558 EFPD, Table 1).
[NotesCreep] "Notes on Creep" (course slides, 2025 revised), slide
             "Warning 2 - Secondary or steady-state stage".
[own]        Derivation written for this module.

Not implemented yet
-------------------
- Stress-enhanced swelling, Ṡ = Ṡ0·(1 + P·σ_hyd) [NotesCreep].
- The (1 + γ1·σ_eq) factor of [Acosta2019] Eq. 9 in the irradiation creep.
- T-dependent incubation and transition doses, [Acosta2019] Eq. 13-14.

Validity of the module as a whole
---------------------------------
- Small strains, isotropic swelling, von Mises (isochoric) creep flow.
- Material parameters live in the card, with their own source and range;
  the functions below only fix the *form* of the laws.
- Written for the 15-15Ti wrapper/cladding of the U_hex_shell work; the T91
  card will reuse the same forms with its own parameters.
"""

import ufl


def _temperature(T, material):
    """Live temperature field, or the card's T_initial in mechanical-only runs.

    Source: z3st convention (same fallback as CreepModel._creep_A).
    Validity: any run; with thermal off the material is isothermal at T_initial.
    """
    if T is None:
        return float(material.get("T_initial", material.get("T_ref", 293.15)))
    return T


def _softplus(x):
    """ln(1 + e^x) in the overflow-safe form max(x, 0) + ln(1 + e^-|x|).

    Source: [own] numerical identity.
    Validity: any real x; used by :func:`volumetric_swelling` with D_tr > 0.
    """
    return ufl.max_value(x, 0.0) + ufl.ln(1.0 + ufl.exp(-abs(x)))


# Void swelling 


def steady_swelling_rate(T, material):
    """g0(T) (1/dpa): steady-state volumetric swelling rate, Gaussian bell in T.

        g0(T) = g0_max · exp(-((T - T_peak) / T_width)²)

    Source: form of [Acosta2019] Eq. 12 (β1·exp(-((T-β2)/β3)²), times Ḋ there,
    per dpa here). [Ye2020] Eq. 2 is the same bell:
    0.3 %/dpa · exp(-4.3e-5 (T - 460 °C)²), i.e. g0_max = 3e-3 1/dpa,
    T_peak = 733.15 K, T_width = 1/sqrt(4.3e-5) = 152.5 K.
    Validity: the T range of the card's calibration. [Acosta2019] gives no
    numerical range nor β values; [Ye2020] states no range and applies it at
    560 °C. Away from the calibration temperature the bell is an extrapolation.
    """
    g0_max = float(material["swelling_g0_max"])
    T_peak = float(material["swelling_T_peak"])
    T_width = float(material["swelling_T_width"])
    return g0_max * ufl.exp(-(((T - T_peak) / T_width) ** 2))


def swelling_rate_per_dpa(T, D, material):
    """dg/dD (1/dpa): logistic onset of g0(T) after the incubation dose.

        dg/dD = g0(T) / (1 + exp((D_inc - D) / D_tr))        (D_tr > 0)
        dg/dD = g0(T) · H(D - D_inc)                         (D_tr = 0)

    Source: [Acosta2019] Eq. 11, [Acosta2021] Eq. 1 (logistic form, symbols
    Δ_g = D_inc, δ_g = D_tr). D_tr = 0 is the bilinear law of [Ye2020] Eq. 2.
    The tanh form 1/(1 + e^-x) = (1 + tanh(x/2))/2 is [own], to avoid overflow
    for a sharp onset.
    Validity: [Acosta2019] uses it for the austenitic cladding of French SFRs
    (DOMAJEUR2), with T-dependent Δ_g and δ_g that are constants here. Used by
    :func:`irradiation_creep` only at dt = 0.
    """
    T = _temperature(T, material)
    g0 = steady_swelling_rate(T, material)
    D_inc = float(material["swelling_D_inc"])
    D_tr = float(material.get("swelling_D_tr", 0.0))
    if D_tr > 0.0:
        return g0 * 0.5 * (1.0 + ufl.tanh((D - D_inc) / (2.0 * D_tr)))
    return ufl.conditional(ufl.gt(D, D_inc), g0, 0.0)


def volumetric_swelling(T, D, material):
    """g = ΔV/V (-) at dose D (dpa), closed form at constant T.

        g(D) = g0 · D_tr · [softplus((D - D_inc)/D_tr) - softplus(-D_inc/D_tr)]
        g(D) = g0 · max(D - D_inc, 0)                         (D_tr = 0)

    Source: the D_tr > 0 expression is [own], the exact integral in D of the
    logistic swelling rate dg/dD = g0 / (1 + exp((D_inc - D)/D_tr)) of
    [Acosta2019] Eq. 11 / [Acosta2021] Eq. 1 (symbols Δ_g = D_inc,
    δ_g = D_tr) at constant T; the integral is not printed in the paper.
    The D_tr = 0 limit is [Ye2020] Eq. 2:
    ΔV/V (%) = 0.3·exp[-4.3e-5 (T - 460)²]·(d - 50) for d > 50 dpa, 0 below.
    Validity: evaluated from the *current* T and D, so path-independent: exact
    only where T stays constant in time at the material point. D_inc and D_tr
    are constants ([Acosta2019] makes them T-dependent). [Ye2020]: 316Ti CW,
    case at 560 °C up to 75 dpa; beyond that it is an extrapolation.
    """
    T = _temperature(T, material)
    g0 = steady_swelling_rate(T, material)
    D_inc = float(material["swelling_D_inc"])
    D_tr = float(material.get("swelling_D_tr", 0.0))
    if D_tr > 0.0:
        return g0 * D_tr * (
            _softplus((D - D_inc) / D_tr) - _softplus(-D_inc / D_tr)
        )
    return g0 * ufl.max_value(D - D_inc, 0.0)


def void_swelling(T, material, model=None, dim=3):
    """Isotropic void-swelling eigenstrain ε* = ((1 + g)^(1/3) - 1) I.

    Card hook: ``eigenstrain: materials.steel_irradiation.void_swelling``,
    added by MechanicalModel.eigenstrain to the thermal ε*.

    Source: (1 + ε_swe)³ = 1 + g, [Acosta2019] Eq. 15, [Acosta2021] Eq. 2;
    isotropy of void swelling [Acosta2021] §2.1. [Ye2020] Eq. 8 uses the
    linearised ε = g/3.
    Validity: isotropic swelling, small strains, no stress dependence; g from
    :func:`volumetric_swelling` (constant-T closed form). Zero until
    ``model.dose`` exists (no ``dpa_rate`` on any card).
    Verified: cases/verification/fuel_assembly/steel_swelling (swelling only) and
    steel_irradiation (coupled with creep).
    """
    I = ufl.Identity(dim)
    D = getattr(model, "dose", None)
    if D is None:
        return 0.0 * I
    g = volumetric_swelling(T, D, material)
    return ((1.0 + g) ** (1.0 / 3.0) - 1.0) * I


# --. irradiation creep --..


def irradiation_creep_uncoupled(T, material, model=None):
    """Simplified irradiation-creep compliance C = B0·Ḋ (1/(Pa·s)), used by
    CreepModel as ε̇_irr = C·σ_eq (card hook
    ``creep_irr: materials.steel_irradiation.irradiation_creep_uncoupled``).

    It is :func:`irradiation_creep` with D_sc = 0: the creep rate does not
    depend on swelling, so swelling and creep can be studied separately.

    Source: the zero-swelling term B0·σ̄·φ of the steady-state law of
    [NotesCreep] (B0 "zero-swelling creep rate coefficient"), written per dpa
    as the first term of [Ye2020] Eq. 3, Δε_eq/σ_eq = 1.1e-6·Δd (σ in MPa).
    Validity: steady-state irradiation creep, linear in σ, isochoric; no T
    dependence (T is not used). Neglecting the swelling-driven term is accurate
    only while swelling is negligible, i.e. below the incubation dose; above it
    the term D_sc·ġ dominates for austenitic steels (e.g. [Ye2020] values).
    [Ye2020]: 316Ti CW, ≤ 75 dpa, 560 °C.
    Verified: cases/verification/fuel_assembly/steel_irradiation_creep.
    """
    Ddot = getattr(model, "dose_rate", None)
    if Ddot is None:
        return 0.0
    B0 = float(material.get("creep_irr_B0", 0.0))
    return B0 * Ddot


def irradiation_creep(T, material, model=None):
    """Irradiation-creep compliance with the swelling-driven term (1/(Pa·s)),
    used by CreepModel as ε̇_irr = C·σ_eq (card hook
    ``creep_irr: materials.steel_irradiation.irradiation_creep``).

        ε̇_eq = (B0 · Ḋ + D_sc · ġ) · σ_eq
        C    = B0·Ḋ + D_sc·[g(D) - g(D - Ḋ·dt)] / dt          (dt > 0)
        C    = Ḋ·(B0 + D_sc·dg/dD)                             (dt = 0)

    D is the end-of-step dose (spine.update_state runs before the solve), so
    g(D - Ḋ·dt) is the start-of-step value.

    Source: the law ε̇_irrad = B0·σ̄·φ + D·Ṡ·σ̄ of [NotesCreep] (B0 zero-swelling
    creep coefficient, D swelling-enhanced creep coefficient, Ṡ swelling rate),
    written per dpa as in [Ye2020] Eq. 3: Δε_eq/σ_eq = 1.1e-6·Δd + 1.2e-2·Δε^S
    (σ in MPa). [Acosta2019] Eq. 9 has the same two terms plus a (1 + γ1 σ)
    factor, not included here. The step-averaged ġ is [own]: it makes the
    backward-Euler creep increment exact at constant T and σ, also when the
    incubation dose falls inside a step (a pointwise rate switches on for the
    whole step at Gauss points where D exceeds D_inc by round-off).
    Validity: steady-state irradiation creep, linear in σ, isochoric; g from
    :func:`volumetric_swelling` (constant-T closed form). [Ye2020]: 316Ti CW,
    ≤ 75 dpa, 560 °C; Ye does not state whether Δε^S is volumetric or linear,
    while this function couples to the volumetric g.
    Verified: cases/verification/fuel_assembly/steel_irradiation (coupled case).
    """
    Ddot = getattr(model, "dose_rate", None)
    D = getattr(model, "dose", None)
    if Ddot is None or D is None:
        return 0.0
    B0 = float(material.get("creep_irr_B0", 0.0))
    D_sc = float(material.get("creep_irr_Dsc", 0.0))
    dt = float(getattr(model, "dt", 0.0) or 0.0)
    if dt > 0.0:
        g_new = volumetric_swelling(T, D, material)
        g_old = volumetric_swelling(T, D - Ddot * dt, material)
        return B0 * Ddot + D_sc * (g_new - g_old) / dt
    return Ddot * (B0 + D_sc * swelling_rate_per_dpa(T, D, material))
