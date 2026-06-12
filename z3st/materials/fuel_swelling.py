# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: burnup-driven swelling eigenstrains for fissile materials
# Author: Giovanni Zullo
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Fuel swelling as a state-dependent eigenstrain (the eigenstrain bus).

A fuel material card names one of these via ``eigenstrain:
materials.fuel_swelling.<name>``; ``spine.load_materials`` resolves it to
``_eigenstrain_func`` and ``MechanicalModel.eigenstrain`` adds its return to the
total inelastic strain ε*, so the eigenstress −ℂ:ε* needs no further change and
the Newton tangent stays automatic.

Signature (the eigenstrain-bus contract)::

    fn(T, material, model=None, dim=3) -> UFL tensor (dim x dim)

- ``T``        : temperature field (UFL) or None.
- ``material`` : the material card dict (read parameters from it).
- ``model``    : the solver/spine — carries the accumulated burnup field
                 ``model.burnup`` (MWd/kgU), written by ``spine.update_state``.
- ``dim``      : eigenstrain tensor dimension (3 for axisymmetric/2d/3d).

These are the parametric first cut of the "fuel is a material" swelling
behaviour; a mechanistic solid+gaseous swelling law (and densification) drops in
behind this same hook later, reading the same burnup/temperature state.
"""

import ufl


def solid_swelling(T, material, model=None, dim=3):
    """Solid fission-product swelling as a burnup-driven volumetric eigenstrain.

        ΔV/V = swelling_rate · bu        (isotropic)
        ε*   = (ΔV/V) / 3 · I            (linear eigenstrain per direction)

    ``bu`` is the accumulated local burnup field (``model.burnup``, MWd/kgU);
    ``swelling_rate`` [ (ΔV/V) per MWd/kgU ] comes from the material card. The
    default 7.0e-4 is the solid fraction of the textbook ~1 %% ΔV/V per
    10 GWd/tU total-swelling rule of thumb (note GWd/tU ≡ MWd/kgU); replace with
    a mechanistic correlation behind this same hook.

    AD-transparent: the burnup field is a fixed coefficient (not the displacement
    unknown), so the eigenstress is a source term and the u-tangent is unchanged.
    """
    I = ufl.Identity(dim)
    rate = float(material.get("swelling_rate", 7.0e-4))   # (ΔV/V) per MWd/kgU
    bu = getattr(model, "burnup", None)
    if bu is None:
        return 0.0 * I
    return (rate * bu / 3.0) * I


def solid_gas_densification(T, material, model=None, dim=3):
    """Combined solid + gaseous swelling and early-life densification.

        ΔV/V = rate_s · bu                                        (solid FP)
             + rate_g · bu · S(T)                                 (gaseous FP)
             − d0 · (1 − exp(−bu / bu_d))                         (densification)
        ε*   = (ΔV/V) / 3 · I

    - Solid swelling: as :func:`solid_swelling` (card ``swelling_rate``).
    - Gaseous swelling: bubble swelling activates with temperature through a
      smooth sigmoid S(T) = 1/(1 + exp(−(T − T_on)/w)) — negligible in cold
      rim, growing toward the pellet centre. Cards ``gas_swelling_rate``
      (default 4.0e-4 ΔV/V per MWd/kgU at full activation), ``gas_T_onset``
      (default 1200 K), ``gas_T_width`` (default 150 K).
    - Densification: in-pile sintering removes as-fabricated porosity early in
      life, ΔV/V → −d0 with a burnup constant bu_d. Cards ``densification_dv``
      (default 0.010, i.e. 1 % ΔV/V recovered) and ``densification_bu``
      (default 2.0 MWd/kgU — ~95 % complete by 6 MWd/kgU). This re-opens the
      gap slightly before swelling closes it — the classic gap-closure shape.

    All terms are UFL expressions in the burnup field (and T for the gaseous
    sigmoid); both are fixed coefficients w.r.t. the displacement unknown, so
    the eigenstress stays a source term and the u-tangent is unchanged.
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
