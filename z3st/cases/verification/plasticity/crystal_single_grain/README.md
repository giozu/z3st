# Single-crystal viscoplasticity with automatic differentiation

## Overview

Single-crystal viscoplasticity on one slip system, implemented as a custom
constitutive law (`plasticity: mode: custom`). The Newton Jacobian is obtained by
UFL symbolic differentiation. The result is checked against the semi-analytical
saturation stress.

In the local suite (carries a gold) and in `cases_ci.txt`.

- FCC slip system (111)[0-11], Schmid factor μ = 0.408
- Power-law viscoplasticity: γ̇ = γ₀ |τ/g₀|ⁿ sign(τ)
- Backward Euler integration of the plastic strain

## Physical model

### Slip system

```
Slip plane:     {111}
Slip direction: <110>
Active system:  (111)[0-11]
Schmid factor:  μ = 0.408248 (for z-axis loading)
```

### Constitutive equations

```
σ = C : (ε_total - ε_p)                    [Stress-strain relation]
ε̇_p = γ̇ · P                                [Plastic strain rate]
γ̇ = γ₀ |τ/g₀|ⁿ sign(τ)                     [Power law slip rate]
τ = σ : P                                   [Resolved shear stress]
P = ½(m⊗n + n⊗m)                            [Schmid tensor]
```

### Time integration (backward Euler)

```
ε_p^{n+1} = ε_p^n + Δt · ε̇_p^{n+1}
```

γ̇ is evaluated at the current stress, so the update is implicit and solved by Newton's method.

### Material parameters (`single_crystal.yaml`, case-local)

| Parameter | Value | Description |
|-----------|-------|-------------|
| E | 200 GPa | Young's modulus |
| ν | 0.3 | Poisson's ratio |
| g₀ | 200 MPa | Slip resistance (CRSS) |
| γ₀ | 0.001 s⁻¹ | Reference slip rate |
| n | 5 | Power law exponent |

The card points to `single_crystal_law.single_crystal_stress` (case-local
`single_crystal_law.py`).

## Geometry and loading

Mesh:
- 1×1×1 m cube
- Transfinite hexahedral mesh: 64 elements, 125 nodes (Gmsh)

Boundary conditions (`boundary_conditions.yaml`):
- Bottom face (z=0): fixed displacement [0, 0, 0]
- Top face (z=1): prescribed displacement in z
- Loading: 41 steps over 2 s, from 0 to 1 % strain (ε̇ = 0.005 s⁻¹)

Solver: Newton (SNES), `direct_mumps`, `rtol: 1e-7`. The stored log shows 1 to 2
staggered iterations per step.

## Semi-analytical solution

Under constant strain rate the stress saturates when the plastic strain rate
equals the total strain rate.

At steady state (dσ/dt = 0):
```
ε̇_total = ε̇_plastic
ε̇_total = m · γ̇ = m · γ₀ · (m·σ/g₀)ⁿ
```

Solving for σ:
```
σ_sat = (g₀/m) · (ε̇_total / (m·γ₀))^(1/n)
```

With the case parameters:
```
ε̇_total = 0.005 s⁻¹
m = 0.408248
g₀ = 200 MPa
γ₀ = 0.001 s⁻¹
n = 5

σ_sat = (200/0.408) × (0.005/(0.408×0.001))^(1/5)
σ_sat = 808.6 MPa
```

Z3ST result (gold): σ_zz = 780.9 MPa at ε_zz = 0.01, 3.4 % below σ_sat.
The resolved shear stress reaches g₀ at σ_zz = g₀/m ≈ 490 MPa.

## Automatic differentiation

With n = 5 the slip-rate derivative is
```
dγ̇/dτ = (n·γ₀/g₀) · (τ/g₀)^(n-1)
```
At τ = 2·g₀ it is 16 times its value at τ = g₀.

The Jacobian is obtained from UFL:
```python
# In single_crystal_law.py
tau_var = ufl.variable(tau)  # Mark as differentiation variable
gamma_dot = gamma0 * (abs(tau_var/g0))**n_pow * ufl.sign(tau_var)

# UFL computes ∂γ̇/∂τ symbolically when assembling the Jacobian
```

## Running the case

```bash
./Allrun
```

Step by step:
```bash
./Allclean                        # clean previous results
python3 -m z3st                   # run
python3 non-regression.py         # checks and stress-strain plot
python3 visualize_material_law.py # material law plots
```

## Schmid factor printout

At start-up `single_crystal_law.py` prints:
```
======================================================================
CRYSTAL PLASTICITY - SCHMID FACTOR CALCULATION
======================================================================
Slip system: (111)[0-11]
  Plane normal n = [0.577, 0.577, 0.577]
  Slip direction m = [0.000, 0.707, -0.707]

Loading direction: e_z = [0, 0, 1]

Schmid factor calculations:
  Method 1 (direct):        μ = |m·e_z| × |n·e_z| = 0.408248
  Method 2 (tensor P_zz):   μ = |P_zz|           = 0.408248
  Method 3 (full tensor):   μ = |P:σ_zz|         = 0.408248

For uniaxial stress σ_zz:
  Resolved shear stress: τ = 0.408248 × σ_zz
======================================================================
```

## Non-regression checks

`non-regression.py` writes `output/stress_strain_curve.png` and checks, with
analytic tolerance 0.25:

| Metric | Reference | Gold value | Relative error |
|--------|-----------|------------|----------------|
| `sigma_zz_final` | σ_sat = 808.6 MPa | 780.9 MPa | 3.4 % |
| `epsilon_zz_final` | 0.01 | 0.01 | 7.6e-13 |
| `saturation_convergence` | 0 % | 3.42 % | 3.4 % |

It also compares every metric with `output/non-regression_gold.json`.

![Stress-Strain Curve](output/stress_strain_curve.png)

## Files

| File | Purpose |
|------|---------|
| `single_crystal_law.py` | Crystal plasticity constitutive model (case-local) |
| `single_crystal.yaml` | Material parameters (case-local) |
| [`plasticity_model.py`](../../../../models/plasticity_model.py) | History variable management |
| `input.yaml` | Simulation configuration |
| `non-regression.py` | Verification and plotting |
| `visualize_material_law.py` | γ̇ against τ, its derivative, n=5 against n=1; writes `output/material_law_visualization.png` |

### History variables

`PlasticityModel` holds:
- `ep_n`: plastic strain tensor at the previous converged step
- `ep`: plastic strain tensor at the current step
- `p_n`: cumulative plastic strain (scalar) at the previous step
- `p`: cumulative plastic strain at the current step

Update sequence:
1. Newton converges and updates `u`.
2. `get_cp_internal_variables()` computes `ep_new = ep_old + Δt·ε̇_p`.
3. `PlasticityModel` sets `ep_n ← ep_new`.
4. The next step uses the updated `ep_n`.

## References

1. Asaro, R. J., & Rice, J. R. (1977). Strain localization in ductile single crystals. *Journal of the Mechanics and Physics of Solids*, 25(5), 309-338.
2. Perzyna, P. (1966). Fundamental problems in viscoplasticity. *Advances in Applied Mechanics*, 9, 243-377.
3. Logg, A., Mardal, K. A., & Wells, G. (2012). *Automated solution of differential equations by the finite element method*. Springer.
4. Baratta, I. A., et al. (2023). DOLFINx: The next generation FEniCS problem solving environment. *Zenodo*.

## Authors

- Giovanni Zullo

Verified against the blessed gold with Z3ST 0.4.0 and FEniCSx 0.11.0.
