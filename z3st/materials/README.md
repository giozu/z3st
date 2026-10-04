# Material cards and property functions

The cards in this directory hold representative values chosen for the
demonstration and verification cases distributed with Z3ST. They are not
qualified design data and are not traceable to a materials handbook unless a
source is named below. Check every value against your own data before using a
card outside these cases. Many cases carry their own card in the case
directory instead of using one from here.

## How a card works

A card is a YAML file of properties in SI units (exception: `oxide.yaml`, see
below). A property is either a number or the dotted path of a Python function:

```yaml
E: 170e+9                 # a number
k: materials.ceramic.k    # a function in materials/ceramic.py
```

When the case is loaded, `Spine.resolve_function` imports the module with
`importlib` and fetches the function. The function receives the temperature
field (and, if its signature asks for them, the material card and the model)
and returns a UFL expression that enters the weak form. Paths resolve against
the `z3st/` package directory, which `z3st/__main__.py` puts on `sys.path`,
and against the working directory, which `python -m z3st` puts there, so a
module in the case directory can be named the same way (as
`verification/plasticity/crystal_single_grain` does with its stress function).
Keys resolved this way: `E`, `nu`, `k`, `Gc`, `eigenstrain`,
`radial_profile`, `axial_profile`, and `stress_function` (resolved by the plasticity
model for `constitutive: custom`). `k` can also be a block that selects a
conductivity model (`type: neural_network`, `gpr` or `magni`).

## Cards

| Card | Contents |
|---|---|
| `steel.yaml` | generic steel, E = 200 GPa, `yield_strength` 200 MPa, `hardening_modulus` 10 GPa for J2 |
| `vessel_steel.yaml` | vessel steel with gamma heating 2 MW/m³ and attenuation 24 1/m |
| `vessel_steel_0.yaml` | the same steel without gamma heating |
| `high_carbon_steel.yaml` | E = 210 GPa, Gc = 2700 J/m² (the values of the SENS/SENT benchmarks of Ambati et al.) |
| `austenitic_steel.yaml` | austenitic stainless steel (AISI 304/316 class) |
| `martensitic_steel.yaml` | martensitic stainless steel (AISI 410/420, 12Cr class) |
| `T91.yaml` | T91 steel, constant properties |
| `15_15Ti.yaml` | 15-15Ti steel; the linear k(T) note in the card cites "Homework 2024-2025, NDT" and is not used by the code |
| `zircaloy.yaml` | Zircaloy-4, E = 99.3 GPa, k = 17 W/m K |
| `uo2.yaml` | UO₂ for the quenched-pellet case: E = 205 GPa, k = 5 W/m K, σc = 1 GPa, initial temperature 1023.15 K |
| `mox_magni.yaml` | MA-MOX with the conductivity of Magni et al. (`magni_mox_thermal.k`) and its composition parameters |
| `ceramic.yaml` | fissile ceramic with `k: materials.ceramic.k` |
| `oxide.yaml` | oxide in micron-based units (pW, pJ, µm), with `k` and a grain-boundary-weakened `Gc(x)` from `oxide.py` |
| `plastic.yaml` | HDPE-like polymer; its `sigma_y` and `plastic` keys are not read by the plasticity model, and the commented `constitutive: voigt` block has no implementation |
| `lead.yaml`, `h2o.yaml` | coolant properties (cp, ρ, viscosity, k) |

## Python property modules

| Module | Functions |
|---|---|
| `ceramic.py` | `k(T)`, a constant 2.5 W/m K written as a function |
| `oxide.py` | `k(T)` (constant, micron units), `Gc(mesh)` and `Gc_numpy(y)`, a tanh transition from a grain-boundary to a bulk toughness |
| `fuel_thermal.py` | `k(T)`, UO₂ modified NFI correlation (Ohira and Itagaki 1997, as adopted in FRAPCON-3) at zero burnup, 95 % TD; burnup degradation not included |
| `magni_mox_thermal.py` | `k(T, ...)` and numpy versions, MA-MOX conductivity of Magni et al. (INSPYRE coefficients) |
| `zircaloy_E.py` | `E(T)`, currently a constant 99.3 GPa written as a UFL expression |
| `fuel_swelling.py` | `solid_gas_densification`, solid and gaseous swelling with early-life densification as an eigenstrain driven by burnup |
| `sciantix_swelling.py` | eigenstrains from the SCIANTIX gaseous-swelling field, alone or with solid swelling and densification |
| `fuel_profiles.py` | radial and axial power form factors: `rim_peaking`, `chopped_cosine`, `tabulated_axial`, `olander_plutonium_redistribution` |
