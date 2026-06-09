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
