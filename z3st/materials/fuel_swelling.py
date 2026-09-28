# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.4.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Fuel swelling as a state-dependent eigenstrain.

A fuel material card names one of these via ``eigenstrain:
materials.fuel_swelling.<name>``; ``spine.load_materials`` resolves it to
``_eigenstrain_func`` and ``MechanicalModel.eigenstrain`` adds its return to the
total inelastic strain ε*.

Signature::

    fn(T, material, model=None, dim=3) -> UFL tensor (dim x dim)

- ``T``        : temperature field (UFL) or None.
- ``material`` : the material card dict (read parameters from it).
- ``model``    : the solver/spine — carries the accumulated burnup field
                 ``model.burnup`` (MWd/kgU), written by ``spine.update_state``.
- ``dim``      : eigenstrain tensor dimension (3 for axisymmetric/2d/3d).

"""

import ufl

def solid_gas_densification(T, material, model=None, dim=3):
    """Combined solid + gaseous swelling and early-life densification.

        ΔV/V = rate_s · bu                                        (solid FP)
             + rate_g · bu · S(T)                                 (gaseous FP)
             - d0 · (1 - exp(- bu / bu_d))                        (densification)
        ε*   = (ΔV/V) / 3 · I

    - Solid swelling: as :func:`solid_swelling` (card ``swelling_rate``).
    - Gaseous swelling: activates with temperature through a smooth sigmoid 
      S(T) = 1/(1 + exp(−(T − T_on)/w)).
      Cards ``gas_swelling_rate``, ``gas_T_onset``, ``gas_T_width``.
    - Densification: Cards ``densification_dv``and ``densification_bu``. 

    All terms are UFL expressions in the burnup field (and T); both are fixed
    coefficients w.r.t. the displacement unknown.
    """
    I = ufl.Identity(dim)
    bu = getattr(model, "burnup", None)
    if bu is None:
        return 0.0 * I

    rate_s = float(material.get("swelling_rate", 7.0e-4))
    rate_g = float(material.get("gas_swelling_rate", 4.0e-4))
    T_on = float(material.get("gas_T_onset", 1200.0))
    width = float(material.get("gas_T_width", 150.0))
    d0 = float(material.get("densification_dv", 0.010))
    bu_d = float(material.get("densification_bu", 2.0))

    dv_solid = rate_s * bu
    if T is None:
        dv_gas = 0.0 * bu
    else:
        S = 1.0 / (1.0 + ufl.exp(-(T - T_on) / width))
        dv_gas = rate_g * bu * S
    dv_dens = -d0 * (1.0 - ufl.exp(-bu / bu_d))

    return ((dv_solid + dv_gas + dv_dens) / 3.0) * I


def uzrh_fission_product_swelling(T, material, model=None, dim=3):
    """Fission-product swelling of U-ZrH (TRIGA) fuel.

        ΔV/V = rate · bu
        ε*   = (ΔV/V) / 3 · I

    Card ``swelling_rate`` in 1/(MWd/kgU), default 2.0e-3 from Simnad &
    Konings, Comprehensive Nuclear Materials vol. 3, §3.12.3.3: hydride fuel
    swells at 3 % per %FIMA, which the same section converts to ΔV/V = 0.2 %
    per MWd/kg equivalent-oxide burnup — three times the 0.07 % of UO2.

    One rate covers all of it. The 3 %/%FIMA is the *total* measured swelling:
    solid fission products, agglomeration of fission gases, and the saturable
    nucleation of irradiation vacancies into voids. Splitting it into a solid
    and a gaseous term the way :func:`solid_gas_densification` does for oxide
    would double-count, and the gaseous sigmoid would be wrong here anyway:
    below ~750 °C the volatiles stay in the hydride (§3.12.3.4), and above it
    the release mechanisms are not those of UO2.

    Two effects are deliberately absent:

    - the offset swelling that precedes the constant-rate regime — the rate is
      the slope *after* it, so this overestimates at very low burnup;
    - the hydride expansion (ΔL/L)_H/Zr = 0.027 (H/Zr − 1.6) of Table 3, which
      needs a spatially varying H/Zr, i.e. a hydrogen field. With a uniform
      H/Zr it is a fabrication strain, not an irradiation eigenstrain. It
      belongs here once a hydrogen-redistribution model provides the field.
    """
    I = ufl.Identity(dim)
    bu = getattr(model, "burnup", None)
    if bu is None:
        return 0.0 * I

    rate = float(material.get("swelling_rate", 2.0e-3))
    return (rate * bu / 3.0) * I
