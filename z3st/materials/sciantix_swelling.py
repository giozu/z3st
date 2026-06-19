# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
SCIANTIX gaseous swelling as a state-dependent eigenstrain.

The mechanical side of the Z3ST <-> SCIANTIX coupling. SCIANTIX runs pointwise
(one integration point per dof, see ``z3st.coupling.sciantix``) and its total
gaseous swelling ΔV/V is written each step into the dolfinx Function
``model.gas_swelling`` by ``spine.update_state``. This callable reads that field
and returns the isotropic eigenstrain::

    ε* = (ΔV/V) / 3 · I

like ``fuel_swelling.solid_swelling`` reads ``model.burnup`` — so the
eigenstress -C:ε* needs no further change and the Newton tangent stays automatic.

A fuel material card opts in with::

    eigenstrain: materials.sciantix_swelling.gaseous_swelling

together with ``models.fission_gas.enabled: true`` (which builds the field and
drives SCIANTIX). To combine SCIANTIX gaseous swelling with a solid-FP swelling
term, use :func:`with_solid` (adds the burnup-linear solid law on top of the
SCIANTIX gaseous field).

Signature::

    fn(T, material, model=None, dim=3) -> UFL tensor (dim x dim)
"""

import ufl


def gaseous_swelling(T, material, model=None, dim=3):
    """Isotropic eigenstrain from the SCIANTIX gaseous-swelling field.

    Reads ``model.gas_swelling`` (a dolfinx Function on the thermal space, total
    ΔV/V from intra- + intergranular gas). Returns ``0`` if the coupling is off
    (no field), so a card may carry this callable harmlessly when SCIANTIX is
    disabled.
    """
    del T, material   # bus-contract args, unused here (depends only on the field)
    I = ufl.Identity(dim)
    gs = getattr(model, "gas_swelling", None)
    if gs is None:
        return 0.0 * I
    return (gs / 3.0) * I


def with_solid(T, material, model=None, dim=3):
    """SCIANTIX gaseous swelling plus the burnup-linear solid-FP swelling.

    ε* = [ (ΔV/V)_gas(SCIANTIX) + rate_s · bu ] / 3 · I

    The gaseous term is the SCIANTIX field; the solid term is the standard
    burnup-linear law (card ``swelling_rate`` [1/(MWd/kgU)], reads ``model.burnup``).
    """
    eps = gaseous_swelling(T, material, model=model, dim=dim)
    bu = getattr(model, "burnup", None)
    if bu is not None:
        rate_s = float(material.get("swelling_rate", 7.0e-4))
        eps = eps + (rate_s * bu / 3.0) * ufl.Identity(dim)
    return eps


def with_solid_densification(T, material, model=None, dim=3):
    """SCIANTIX gaseous swelling plus burnup-linear solid-FP swelling and early-life
    densification — the SCIANTIX equivalent of ``fuel_swelling.solid_gas_densification``.

        ΔV/V = (ΔV/V)_gas(SCIANTIX)                            (intra+intergranular)
             + rate_s · bu                                     (solid FP)
             - d0 · (1 - exp(- bu / bu_d))                     (densification)
        ε*   = (ΔV/V) / 3 · I

    Use this in place of ``solid_gas_densification`` to take the gaseous term from the
    SCIANTIX field instead of the analytic sigmoid; the solid and
    densification terms (cards ``swelling_rate``, ``densification_dv``,
    ``densification_bu``) are unchanged and stay UFL expressions in ``model.burnup``,
    so the u-tangent is untouched. SCIANTIX runs with ``iDensification = 0`` — Z3ST
    owns densification — so the two do not double-count.
    """
    eps = gaseous_swelling(T, material, model=model, dim=dim)
    bu = getattr(model, "burnup", None)
    if bu is None:
        return eps
    I = ufl.Identity(dim)
    rate_s = float(material.get("swelling_rate", 7.0e-4))
    d0 = float(material.get("densification_dv", 0.010))
    bu_d = float(material.get("densification_bu", 2.0))
    dv_solid = rate_s * bu
    dv_dens = -d0 * (1.0 - ufl.exp(-bu / bu_d))
    return eps + ((dv_solid + dv_dens) / 3.0) * I
