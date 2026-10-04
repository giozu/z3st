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


def uzrh_isotropic_fission_product_eigenstrain(fima, dim=3):
    """Existing U-ZrH correlation, independent of FIMA provenance.

    DeltaV/V = 3*FIMA; infinitesimal isotropic eigenstrain = FIMA*I.
    """
    return fima * ufl.Identity(dim)


def uzrh_fission_product_swelling(T, material, model=None, dim=3):
    """SNAP/Olander literature correlation for U-ZrH fission-product swelling.

    Read ``model.burnup`` in real MWd/kgU and convert to the dimensionless
    fraction FIMA = fissions / (initial U + Zr atoms), excluding hydrogen::

        r = N_U / N_Zr = w_U * (M_Zr + x*M_H) / ((1 - w_U)*M_U)
        FIMA = bu * 8.64e10 * M_U / (N_A * E_f) * r / (1 + r)
        DeltaV/V = 3 * FIMA
        epsilon_sw = FIMA * I

    Required cards: ``heavy_metal_fraction`` is kgU/kg total U-ZrH fuel;
    ``hydrogen_zirconium_ratio`` is the atomic ratio x = H/Zr. Assume a
    mixture of U and ZrH_x, with M_U = 0.238 kg/mol (isotopic approximation).
    ``swelling_energy_per_fission_MeV`` defaults to 200 MeV of deposited
    thermal energy per fission. Density cancels in this conversion.

    The approximately 3% volumetric swelling per %FIMA slope comes from
    U-ZrH/SNAP data discussed by Olander et al., "Uranium-zirconium hydride
    fuel properties" (2009). Its use for H/Zr = 1.0 is an approximation to
    verify against specific irradiation data, not a correlation validated
    specifically for TRIGA FE101.

    This represents the linear contribution after offset swelling; no offset
    or onset threshold is modeled. Hydrogen-redistribution expansion is also
    absent. The former ``swelling_rate = 2.0e-3`` used an equivalent-oxide
    burnup basis and is not directly compatible with real MWd/kgU; this
    function no longer reads that card.
    """
    # Explicit native ownership: never silently fall back to the BU estimate.
    source = material.get("fima_source", "legacy_burnup")
    if source == "native_openmc":
        fima = getattr(model, "fima_native", None)
        if fima is None:
            raise ValueError("native_openmc FIMA field is unavailable")
        return uzrh_isotropic_fission_product_eigenstrain(fima, dim)
    if source != "legacy_burnup":
        raise ValueError(f"Unknown FIMA source: {source}")
    I = ufl.Identity(dim)
    bu = getattr(model, "burnup", None)
    if bu is None:
        return 0.0 * I

    w_u = float(material["heavy_metal_fraction"])
    h_zr = float(material["hydrogen_zirconium_ratio"])
    energy_mev = float(material.get("swelling_energy_per_fission_MeV", 200.0))
    if not 0.0 < w_u < 1.0:
        raise ValueError("U-ZrH swelling requires 0 < heavy_metal_fraction < 1")
    if not h_zr >= 0.0 or not energy_mev > 0.0:
        raise ValueError("Invalid H/Zr ratio or energy per fission")

    M_U = 0.238       # kg/mol
    M_ZR = 0.091224   # kg/mol
    M_H = 0.001008    # kg/mol
    N_A = 6.02214076e23  # atoms/mol
    energy_j = energy_mev * 1.602176634e-13  # J/fission

    u_zr = w_u * (M_ZR + h_zr * M_H) / ((1.0 - w_u) * M_U)
    fima = bu * 8.64e10 * M_U / (N_A * energy_j) * u_zr / (1.0 + u_zr)
    return uzrh_isotropic_fission_product_eigenstrain(fima, dim)
