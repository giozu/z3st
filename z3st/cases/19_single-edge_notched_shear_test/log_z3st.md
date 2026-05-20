Info    : Reading 'mesh.msh'...
Info    : 13 entities
Info    : 83428 nodes
Info    : 167479 elements
Info    : Done reading 'mesh.msh'


***

Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
Author: Giovanni Zullo
Version: 0.1.0 (2025)

***



## Description

Z3ST is an open-source framework for the thermo-mechanical modelling
of materials. Built on FEniCSx, it supports transient simulations,
complex geometries, and user-defined boundary conditions.


### Config initializer

  → Geometry            : geometry.yaml
  → Mesh                : mesh.msh
  → Boundary conditions : boundary_conditions.yaml
  → Time steps          : 796
  → Regime              : 2d
  → Models active       :
      thermal    → OFF
      mechanical → ON
      damage     → ON
      cluster    → OFF
      plasticity → OFF
  → Gap conductance     : None (value = 0.0)



### FiniteElementSetup initializer

Mechanical element order: 1
Thermal function space (V_t): FunctionSpace(<Mesh #0>, Basix element (P, triangle, 1, gll_warped, unset, False, float64, []))
Mechanical function space (V_m): FunctionSpace(<Mesh #0>, blocked element (Basix element (P, triangle, 1, gll_warped, unset, False, float64, []), (2,)))
Scalar function space (V_d): FunctionSpace(<Mesh #0>, Basix element (P, triangle, 1, gll_warped, unset, False, float64, []))
Scalar function space (Q): FunctionSpace(<Mesh #0>, Basix element (P, triangle, 0, gll_warped, unset, True, float64, []))
[Solver] initializer
  Applied relaxation factor:
  → Temperature  : 0.9
  → Displacement : 1.0
  → Damage       : 0.8
  Adaptive relaxation disabled


[MechanicalModel] initializer
[MechanicalModel] options loaded from input.yaml:
  solver              : linear
  linear_solver       : direct_mumps
  rtol                : 1e-05
  stag_tol            : 1e-05
  convergence         : rel_norm
DamageModel initializer
Options loaded from input.yaml:
  type                : AT2
  solver              : linear
  linear_solver       : direct_mumps
  rtol                : 1e-05
  stag_tol            : 1e-07
  convergence         : rel_norm
  lc                  : 4e-06
  hybrid_constraint   : True
  split               : star_convex
  gamma_star          : 5.0
[spine.load_materials]
Material loaded: steel
  → k defined as constant: 45.0
  → Gc defined as constant: 2700.0
  - Material 'steel': sigma_c (AT2) from Gc = 2700.00 J/m2
  → constitutive model: lame
  E               → 210000000000.0 (float)
  G               → 80769230769.23077 (float)
  Gc              → 2700.0 (float)
  T_ref           → 293.15 (float)
  alpha           → 1.1e-05 (float)
  bulk_modulus    → 174999999999.99997 (float)
  constitutive_mode → lame (str)
  cp              → 470.0 (float)
  gamma_heating   → 0.0 (float)
  k               → 45.0 (float)
  lmbda           → 121153846153.84615 (float)
  mu_gamma        → 25.0 (float)
  name            → high_carbon_steel (str)
  nu              → 0.3 (float)
  rho             → 7850.0 (float)
  sigma_c         → 3866548242.61899 (float)
[spine.initialize_fields]

Initializing the displacement field...
  Initial u: min=0.00e+00 m, max=0.00e+00 m, mean=0.00e+00 m

Initializing the damage field...



***


### spine - set_boundary_conditions


***



Loading boundary conditions from 'boundary_conditions.yaml'
  **[INFO]** Constant Dirichlet vector (2D) → [0.0, 0.0]
  **[INFO]** Dirichlet mechanical BC on 'steel' → [0.0, 0.0] at region 'ymin'
  **[INFO]** Step-dependent Dirichlet list (2D), length 796
  **[INFO]** Dirichlet mechanical BC on 'steel' → [0.0, 0.0] at region 'ymax'

Setting damage boundary conditions...
  **[INFO]** Dirichlet damage BC on 'steel' → D = 1.0 at region 'crack'
Computing symbolic re# NOTE: the same plot is reproduced by `plot_energy_balance.py` as a
# standalone script that accepts an `OUTPUT_DIR` positional argument. Use
# that script for post-hoc plotting against a saved backup directory
# (e.g. output_starconvex_g00/) without disturbing the live-run state at
# the case root. The block below is preserved here so a single
# `./Allrun` end-to-end pipeline still produces every diagnostic.sult fields (strain, stress, ...)

**[INFO]** Hot-reload of allow-listed input.yaml parameters is active. Edit input.yaml during the run; changes apply at the next step boundary. Allowed keys: damage.{stag_tol,rtol,hybrid_constraint,gamma_star}, mechanical.{stag_tol,rtol}, thermal.{stag_tol,rtol}, solver_settings.{max_iters,relax_*}.


## Step 01/796: t = 0.00e+00 s | LHR = 0.00e+00 W/m

  → dt=0: solving static step / initial condition
Computing symbolic result fields (strain, stress, ...)



***


### spine - solve


***



Current step = 0 | dt = 0.00e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.000e+00
  |ΔD|_∞ = 8.000e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.000e-01
  |ΔD|_∞ = 1.600e-01

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.000e-02
  |ΔD|_∞ = 3.200e-02

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.000e-03
  |ΔD|_∞ = 6.400e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.600e-03
  |ΔD|_∞ = 1.280e-03

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.200e-04
  |ΔD|_∞ = 2.560e-04

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.400e-05
  |ΔD|_∞ = 5.120e-05

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.280e-05
  |ΔD|_∞ = 1.024e-05

Convergence check


#### Iteration 9/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.560e-06
  |ΔD|_∞ = 2.048e-06

Convergence check


#### Iteration 10/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.120e-07
  |ΔD|_∞ = 4.096e-07

Convergence check


#### Iteration 11/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.024e-07
  |ΔD|_∞ = 8.192e-08

Convergence check


#### Iteration 12/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.048e-08
  |ΔD|_∞ = 1.638e-08

Convergence check

**[SUCCESS]** Staggered solver converged in 12 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 0.0000e+00 J
  → Fracture energy : 1.3579e+00 J
  → Total energy    : 1.3579e+00 J


## Step 02/796: t = 4.53e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 1 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.294e-04
  |ΔD|_∞ = 7.928e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.259e-04
  |ΔD|_∞ = 1.586e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.518e-05
  |ΔD|_∞ = 3.171e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.036e-06
  |ΔD|_∞ = 6.342e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.007e-06
  |ΔD|_∞ = 1.268e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.014e-07
  |ΔD|_∞ = 2.537e-07

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.028e-08
  |ΔD|_∞ = 5.074e-08

Convergence check

**[SUCCESS]** Staggered solver converged in 7 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5969e-02 J
  → Fracture energy : 1.3579e+00 J
  → Total energy    : 1.3839e+00 J


## Step 03/796: t = 9.06e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 2 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.000e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.881e-03
  |ΔD|_∞ = 2.369e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.306e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.761e-04
  |ΔD|_∞ = 4.738e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.522e-05
  |ΔD|_∞ = 9.477e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.210e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.504e-05
  |ΔD|_∞ = 1.895e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.210e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.009e-06
  |ΔD|_∞ = 3.791e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.018e-07
  |ΔD|_∞ = 7.581e-07

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.204e-07
  |ΔD|_∞ = 1.516e-07

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.407e-08
  |ΔD|_∞ = 3.033e-08

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0375e-01 J
  → Fracture energy : 1.3580e+00 J
  → Total energy    : 1.4617e+00 J


## Step 04/796: t = 1.36e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 3 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.334e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.164e-03
  |ΔD|_∞ = 4.063e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.328e-04
  |ΔD|_∞ = 8.125e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.266e-04
  |ΔD|_∞ = 1.625e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.334e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.531e-05
  |ΔD|_∞ = 3.250e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.334e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.063e-06
  |ΔD|_∞ = 6.500e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.013e-06
  |ΔD|_∞ = 1.300e-06

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.025e-07
  |ΔD|_∞ = 2.600e-07

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.050e-08
  |ΔD|_∞ = 5.200e-08

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3295e-01 J
  → Fracture energy : 1.3583e+00 J
  → Total energy    : 1.5913e+00 J


## Step 05/796: t = 1.81e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 4 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.502e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.539e-03
  |ΔD|_∞ = 6.051e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.078e-04
  |ΔD|_∞ = 1.210e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.816e-04
  |ΔD|_∞ = 2.420e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.631e-05
  |ΔD|_∞ = 4.841e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.262e-06
  |ΔD|_∞ = 9.682e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.452e-06
  |ΔD|_∞ = 1.936e-06

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.905e-07
  |ΔD|_∞ = 3.873e-07

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.810e-08
  |ΔD|_∞ = 7.746e-08

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.1286e-01 J
  → Fracture energy : 1.3592e+00 J
  → Total energy    : 1.7721e+00 J


## Step 06/796: t = 2.26e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 5 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.003e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.077e-03
  |ΔD|_∞ = 8.562e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.591e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.215e-03
  |ΔD|_∞ = 1.712e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.591e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.431e-04
  |ΔD|_∞ = 3.425e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.862e-05
  |ΔD|_∞ = 6.849e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.389e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.723e-06
  |ΔD|_∞ = 1.370e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.652e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.945e-06
  |ΔD|_∞ = 2.740e-06

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.660e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.889e-07
  |ΔD|_∞ = 5.480e-07

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  

## Step 07/796: t = 2.72e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 6 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.926e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.634e-03
  |ΔD|_∞ = 5.362e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.003e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.269e-04
  |ΔD|_∞ = 1.072e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.003e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.054e-04
  |ΔD|_∞ = 2.145e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.088e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.108e-05
  |ΔD|_∞ = 4.290e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.885e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.215e-06
  |ΔD|_∞ = 8.580e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.885e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.430e-07
  |ΔD|_∞ = 1.716e-06

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.887e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.686e-07
  |ΔD|_∞ = 3.432e-07

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.182e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.372e-08
  |ΔD|_∞ = 6.864e-08

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.9328e-01 J
  → Fracture energy : 1.3624e+00 J
  → Total energy    : 2.0557e+00 J


## Step 08/796: t = 3.17e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 7 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.741e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.151e-03
  |ΔD|_∞ = 4.535e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.170e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.301e-04
  |ΔD|_∞ = 9.071e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.170e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.603e-05
  |ΔD|_∞ = 1.814e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.033e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.721e-05
  |ΔD|_∞ = 3.628e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.033e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.441e-06
  |ΔD|_∞ = 7.256e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.882e-07
  |ΔD|_∞ = 1.451e-06

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.376e-07
  |ΔD|_∞ = 2.903e-07

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.753e-08
  |ΔD|_∞ = 5.805e-08

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.4638e-01 J
  → Fracture energy : 1.3635e+00 J
  → Total energy    : 2.1099e+00 J


## Step 09/796: t = 3.62e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 8 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.602e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.126e-03
  |ΔD|_∞ = 4.509e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.252e-04
  |ΔD|_∞ = 9.018e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.505e-05
  |ΔD|_∞ = 1.804e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.045e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.701e-05
  |ΔD|_∞ = 3.607e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.023e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.402e-06
  |ΔD|_∞ = 7.214e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.023e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.804e-07
  |ΔD|_∞ = 1.443e-06

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.970e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.361e-07
  |ΔD|_∞ = 2.886e-07

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.861e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.722e-08
  |ΔD|_∞ = 5.771e-08

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.0132e-01 J
  → Fracture energy : 1.3647e+00 J
  → Total energy    : 2.1660e+00 J


## Step 10/796: t = 4.08e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 9 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.480e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.252e-03
  |ΔD|_∞ = 4.838e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.503e-04
  |ΔD|_∞ = 9.677e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.031e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.007e-05
  |ΔD|_∞ = 1.935e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.031e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.801e-05
  |ΔD|_∞ = 3.871e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.276e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.603e-06
  |ΔD|_∞ = 7.741e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.760e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.206e-07
  |ΔD|_∞ = 1.548e-06

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.823e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.441e-07
  |ΔD|_∞ = 3.097e-07

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.882e-08
  |ΔD|_∞ = 6.193e-08

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.5798e-01 J
  → Fracture energy : 1.3661e+00 J
  → Total energy    : 2.2241e+00 J


## Step 11/796: t = 4.53e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 10 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.370e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.462e-03
  |ΔD|_∞ = 5.465e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.234e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.924e-04
  |ΔD|_∞ = 1.093e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.918e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.848e-05
  |ΔD|_∞ = 2.186e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.918e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.970e-05
  |ΔD|_∞ = 4.372e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.939e-06
  |ΔD|_∞ = 8.744e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.879e-07
  |ΔD|_∞ = 1.749e-06

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.576e-07
  |ΔD|_∞ = 3.498e-07

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.151e-08
  |ΔD|_∞ = 6.995e-08

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.1624e-01 J
  → Fracture energy : 1.3678e+00 J
  → Total energy    : 2.2840e+00 J


## Step 12/796: t = 4.98e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 11 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.270e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.746e-03
  |ΔD|_∞ = 6.499e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.491e-04
  |ΔD|_∞ = 1.300e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.630e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.098e-04
  |ΔD|_∞ = 2.600e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.630e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.197e-05
  |ΔD|_∞ = 5.199e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.393e-06
  |ΔD|_∞ = 1.040e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.813e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.786e-07
  |ΔD|_∞ = 2.080e-06

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.029e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.757e-07
  |ΔD|_∞ = 4.159e-07

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.630e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.514e-08
  |ΔD|_∞ = 8.319e-08

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.7594e-01 J
  → Fracture energy : 1.3699e+00 J
  → Total energy    : 2.3458e+00 J


## Step 13/796: t = 5.43e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 12 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.181e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.125e-03
  |ΔD|_∞ = 7.848e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.392e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.249e-04
  |ΔD|_∞ = 1.570e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.392e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.250e-04
  |ΔD|_∞ = 3.139e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.500e-05
  |ΔD|_∞ = 6.279e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.126e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.999e-06
  |ΔD|_∞ = 1.256e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.366e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.999e-07
  |ΔD|_∞ = 2.512e-06

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.367e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.000e-07
  |ΔD|_∞ = 5.023e-07

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.502e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.999e-08
  |ΔD|_∞ = 1.005e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0369e+00 J
  → Fracture energy : 1.3725e+00 J
  → Total energy    : 2.4093e+00 J


## Step 14/796: t = 5.89e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 13 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.104e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.635e-03
  |ΔD|_∞ = 9.357e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.270e-04
  |ΔD|_∞ = 1.871e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.293e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.454e-04
  |ΔD|_∞ = 3.743e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.293e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.908e-05
  |ΔD|_∞ = 7.485e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.816e-06
  |ΔD|_∞ = 1.497e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.163e-06
  |ΔD|_∞ = 2.994e-06

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.326e-07
  |ΔD|_∞ = 5.988e-07

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.010e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.653e-08
  |ΔD|_∞ = 1.198e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0987e+00 J
  → Fracture energy : 1.3758e+00 J
  → Total energy    : 2.4744e+00 J


## Step 15/796: t = 6.34e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 14 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.045e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.313e-03
  |ΔD|_∞ = 1.084e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.626e-04
  |ΔD|_∞ = 2.167e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.725e-04
  |ΔD|_∞ = 4.335e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.450e-05
  |ΔD|_∞ = 8.670e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.979e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.900e-06
  |ΔD|_∞ = 1.734e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.979e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.380e-06
  |ΔD|_∞ = 3.468e-06

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.096e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.760e-07
  |ΔD|_∞ = 6.936e-07

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.096e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.520e-08
  |ΔD|_∞ = 1.387e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1608e+00 J
  → Fracture energy : 1.3802e+00 J
  → Total energy    : 2.5410e+00 J


## Step 16/796: t = 6.79e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 15 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.014e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.149e-03
  |ΔD|_∞ = 1.315e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.749e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.030e-03
  |ΔD|_∞ = 2.630e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.060e-04
  |ΔD|_∞ = 5.261e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.119e-05
  |ΔD|_∞ = 1.052e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.238e-06
  |ΔD|_∞ = 2.104e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.870e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.648e-06
  |ΔD|_∞ = 4.208e-06

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.611e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.295e-07
  |ΔD|_∞ = 8.417e-07

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.938e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.590e-08
  |ΔD|_∞ = 1.683e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2224e+00 J
  → Fracture energy : 1.3862e+00 J
  → Total energy    : 2.6085e+00 J


## Step 17/796: t = 7.25e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 16 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.028e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.095e-03
  |ΔD|_∞ = 1.618e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.219e-03
  |ΔD|_∞ = 3.237e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.438e-04
  |ΔD|_∞ = 6.474e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.876e-05
  |ΔD|_∞ = 1.295e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.752e-06
  |ΔD|_∞ = 2.589e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.950e-06
  |ΔD|_∞ = 5.179e-06

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.298e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.901e-07
  |ΔD|_∞ = 1.036e-06

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.298e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.802e-08
  |ΔD|_∞ = 2.072e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2824e+00 J
  → Fracture energy : 1.3941e+00 J
  → Total energy    : 2.6765e+00 J


## Step 18/796: t = 7.70e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 17 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.109e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.015e-03
  |ΔD|_∞ = 1.885e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.403e-03
  |ΔD|_∞ = 3.770e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.806e-04
  |ΔD|_∞ = 7.539e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.796e-20

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.612e-05
  |ΔD|_∞ = 1.508e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.796e-20

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.122e-05
  |ΔD|_∞ = 3.016e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.245e-06
  |ΔD|_∞ = 6.031e-06

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.490e-07
  |ΔD|_∞ = 1.206e-06

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.979e-08
  |ΔD|_∞ = 2.412e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3399e+00 J
  → Fracture energy : 1.4042e+00 J
  → Total energy    : 2.7442e+00 J


## Step 19/796: t = 8.15e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 18 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.263e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.833e-03
  |ΔD|_∞ = 1.981e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.567e-03
  |ΔD|_∞ = 3.963e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.807e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.133e-04
  |ΔD|_∞ = 7.926e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.807e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.267e-05
  |ΔD|_∞ = 1.585e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.256e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.253e-05
  |ΔD|_∞ = 3.170e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.257e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.507e-06
  |ΔD|_∞ = 6.340e-06

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.013e-07
  |ΔD|_∞ = 1.268e-06

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.003e-07
  |ΔD|_∞ = 2.536e-07

Convergence check


#### Iteration 9/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.005e-08
  |ΔD|_∞ = 5.072e-08

Convergence check

**[SUCCESS]** Staggered solver converged in 9 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3947e+00 J
  → Fracture energy : 1.4164e+00 J
  → Total energy    : 2.8111e+00 J


## Step 20/796: t = 8.60e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 19 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.466e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.657e-03
  |ΔD|_∞ = 2.277e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.543e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.731e-03
  |ΔD|_∞ = 4.553e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.463e-04
  |ΔD|_∞ = 9.107e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.925e-05
  |ΔD|_∞ = 1.821e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.385e-05
  |ΔD|_∞ = 3.643e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.770e-06
  |ΔD|_∞ = 7.285e-06

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.540e-07
  |ΔD|_∞ = 1.457e-06

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.121e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.108e-07
  |ΔD|_∞ = 2.914e-07

Convergence check


#### Iteration 9/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.282e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.216e-08
  |ΔD|_∞ = 5.828e-08

Convergence check

**[SUCCESS]** Staggered solver converged in 9 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4468e+00 J
  → Fracture energy : 1.4301e+00 J
  → Total energy    : 2.8769e+00 J


## Step 21/796: t = 9.06e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 20 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.772e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.198e-03
  |ΔD|_∞ = 2.607e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.840e-03
  |ΔD|_∞ = 5.214e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.043e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.679e-04
  |ΔD|_∞ = 1.043e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.117e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.358e-05
  |ΔD|_∞ = 2.086e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.117e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.472e-05
  |ΔD|_∞ = 4.171e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.943e-06
  |ΔD|_∞ = 8.343e-06

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.887e-07
  |ΔD|_∞ = 1.669e-06

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.050e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.177e-07
  |ΔD|_∞ = 3.337e-07

Convergence check


#### Iteration 9/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.050e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.355e-08
  |ΔD|_∞ = 6.674e-08

Convergence check

**[SUCCESS]** Staggered solver converged in 9 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4954e+00 J
  → Fracture energy : 1.4449e+00 J
  → Total energy    : 2.9403e+00 J


## Step 22/796: t = 9.51e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 21 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.01e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.950e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.738e-03
  |ΔD|_∞ = 2.713e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.01e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.548e-03
  |ΔD|_∞ = 5.425e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.01e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.260e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.095e-04
  |ΔD|_∞ = 1.085e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.01e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.260e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.190e-05
  |ΔD|_∞ = 2.170e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.01e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.238e-05
  |ΔD|_∞ = 4.340e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.01e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.476e-06
  |ΔD|_∞ = 8.680e-06

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.01e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.952e-07
  |ΔD|_∞ = 1.736e-06

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.01e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.905e-08
  |ΔD|_∞ = 3.472e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4735e+00 J
  → Fracture energy : 1.4575e+00 J
  → Total energy    : 2.9310e+00 J


## Step 23/796: t = 9.96e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 22 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.02e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.787e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.417e-03
  |ΔD|_∞ = 2.753e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.02e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.083e-03
  |ΔD|_∞ = 5.506e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.02e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.717e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.167e-04
  |ΔD|_∞ = 1.101e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.02e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.724e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.333e-05
  |ΔD|_∞ = 2.202e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.02e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.592e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.667e-06
  |ΔD|_∞ = 4.405e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.02e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.733e-06
  |ΔD|_∞ = 8.809e-06

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.02e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.467e-07
  |ΔD|_∞ = 1.762e-06

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.02e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.933e-08
  |ΔD|_∞ = 3.524e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4578e+00 J
  → Fracture energy : 1.4664e+00 J
  → Total energy    : 2.9242e+00 J


## Step 24/796: t = 1.04e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 23 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.03e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.849e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.171e-03
  |ΔD|_∞ = 1.983e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.03e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.077e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.342e-04
  |ΔD|_∞ = 3.966e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.03e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.077e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.268e-04
  |ΔD|_∞ = 7.932e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.03e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.445e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.537e-05
  |ΔD|_∞ = 1.586e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.03e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.445e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.074e-06
  |ΔD|_∞ = 3.173e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.03e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.015e-06
  |ΔD|_∞ = 6.346e-06

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.03e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.029e-07
  |ΔD|_∞ = 1.269e-06

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.03e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.122e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.059e-08
  |ΔD|_∞ = 2.538e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4516e+00 J
  → Fracture energy : 1.4712e+00 J
  → Total energy    : 2.9228e+00 J


## Step 25/796: t = 1.09e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 24 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.04e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.577e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.830e-03
  |ΔD|_∞ = 1.038e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.04e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.160e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.659e-04
  |ΔD|_∞ = 2.075e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.04e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.014e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.318e-05
  |ΔD|_∞ = 4.151e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.04e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.014e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.464e-05
  |ΔD|_∞ = 8.301e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.04e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.927e-06
  |ΔD|_∞ = 1.660e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.04e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.854e-07
  |ΔD|_∞ = 3.321e-06

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.04e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.171e-07
  |ΔD|_∞ = 6.641e-07

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.04e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.342e-08
  |ΔD|_∞ = 1.328e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4511e+00 J
  → Fracture energy : 1.4736e+00 J
  → Total energy    : 2.9247e+00 J


## Step 26/796: t = 1.13e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 25 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.05e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.847e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.123e-03
  |ΔD|_∞ = 7.311e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.05e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.245e-04
  |ΔD|_∞ = 1.462e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.05e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.490e-05
  |ΔD|_∞ = 2.925e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.05e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.672e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.980e-06
  |ΔD|_∞ = 5.849e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.05e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.672e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.796e-06
  |ΔD|_∞ = 1.170e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.05e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.592e-07
  |ΔD|_∞ = 2.340e-06

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.05e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.184e-08
  |ΔD|_∞ = 4.679e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 7 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4529e+00 J
  → Fracture energy : 1.4748e+00 J
  → Total energy    : 2.9277e+00 J


## Step 27/796: t = 1.18e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 26 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.06e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.187e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.507e-04
  |ΔD|_∞ = 5.192e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.06e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.501e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.501e-04
  |ΔD|_∞ = 1.038e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.06e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.003e-05
  |ΔD|_∞ = 2.077e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.06e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.005e-06
  |ΔD|_∞ = 4.154e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.06e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.201e-06
  |ΔD|_∞ = 8.307e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.06e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.501e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.402e-07
  |ΔD|_∞ = 1.661e-06

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.06e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.501e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.804e-08
  |ΔD|_∞ = 3.323e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 7 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4555e+00 J
  → Fracture energy : 1.4755e+00 J
  → Total energy    : 2.9311e+00 J


## Step 28/796: t = 1.22e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 27 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.07e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.668e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.292e-04
  |ΔD|_∞ = 3.574e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.07e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.058e-04
  |ΔD|_∞ = 7.147e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.07e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.117e-05
  |ΔD|_∞ = 1.429e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.07e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.579e-20

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.234e-06
  |ΔD|_∞ = 2.859e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.07e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.579e-20

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.467e-07
  |ΔD|_∞ = 5.718e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.07e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.693e-07
  |ΔD|_∞ = 1.144e-06

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.07e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.387e-08
  |ΔD|_∞ = 2.287e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 7 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4585e+00 J
  → Fracture energy : 1.4760e+00 J
  → Total energy    : 2.9345e+00 J


## Step 29/796: t = 1.27e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 28 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.08e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.465e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.869e-04
  |ΔD|_∞ = 2.685e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.08e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.737e-05
  |ΔD|_∞ = 5.370e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.08e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.547e-05
  |ΔD|_∞ = 1.074e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.08e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.095e-06
  |ΔD|_∞ = 2.148e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.08e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.037e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.190e-07
  |ΔD|_∞ = 4.296e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.08e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.037e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.238e-07
  |ΔD|_∞ = 8.593e-07

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.08e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.476e-08
  |ΔD|_∞ = 1.719e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 7 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4617e+00 J
  → Fracture energy : 1.4763e+00 J
  → Total energy    : 2.9381e+00 J


## Step 30/796: t = 1.31e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 29 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.09e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.371e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.988e-04
  |ΔD|_∞ = 2.758e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.09e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.976e-05
  |ΔD|_∞ = 5.517e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.09e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.534e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.195e-05
  |ΔD|_∞ = 1.103e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.09e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.944e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.390e-06
  |ΔD|_∞ = 2.207e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.09e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.572e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.781e-07
  |ΔD|_∞ = 4.413e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.09e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.751e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.561e-08
  |ΔD|_∞ = 8.827e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4650e+00 J
  → Fracture energy : 1.4766e+00 J
  → Total energy    : 2.9416e+00 J


## Step 31/796: t = 1.36e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 30 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.1e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.326e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.400e-04
  |ΔD|_∞ = 2.089e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.1e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.245e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.799e-05
  |ΔD|_∞ = 4.179e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.1e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.374e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.598e-06
  |ΔD|_∞ = 8.358e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.1e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.374e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.920e-06
  |ΔD|_∞ = 1.672e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.1e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.337e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.839e-07
  |ΔD|_∞ = 3.343e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.1e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.336e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.679e-08
  |ΔD|_∞ = 6.686e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4684e+00 J
  → Fracture energy : 1.4768e+00 J
  → Total energy    : 2.9452e+00 J


## Step 32/796: t = 1.40e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 31 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.11e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.300e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.974e-04
  |ΔD|_∞ = 1.323e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.11e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.948e-05
  |ΔD|_∞ = 2.646e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.11e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.378e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.897e-06
  |ΔD|_∞ = 5.292e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.11e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.404e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.579e-06
  |ΔD|_∞ = 1.058e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.11e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.450e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.159e-07
  |ΔD|_∞ = 2.117e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.11e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.317e-08
  |ΔD|_∞ = 4.234e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4718e+00 J
  → Fracture energy : 1.4770e+00 J
  → Total energy    : 2.9488e+00 J


## Step 33/796: t = 1.45e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 32 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.12e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.283e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.713e-04
  |ΔD|_∞ = 1.111e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.12e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.426e-05
  |ΔD|_∞ = 2.222e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.12e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.851e-06
  |ΔD|_∞ = 4.444e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.12e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.370e-06
  |ΔD|_∞ = 8.888e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.12e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.740e-07
  |ΔD|_∞ = 1.778e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.12e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.447e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.481e-08
  |ΔD|_∞ = 3.555e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4753e+00 J
  → Fracture energy : 1.4771e+00 J
  → Total energy    : 2.9524e+00 J


## Step 34/796: t = 1.49e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 33 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.13e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.273e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.524e-04
  |ΔD|_∞ = 9.365e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.13e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.048e-05
  |ΔD|_∞ = 1.873e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.13e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.096e-06
  |ΔD|_∞ = 3.746e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.13e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.219e-06
  |ΔD|_∞ = 7.492e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.13e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.450e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.438e-07
  |ΔD|_∞ = 1.498e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.13e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.450e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.877e-08
  |ΔD|_∞ = 2.997e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4788e+00 J
  → Fracture energy : 1.4772e+00 J
  → Total energy    : 2.9560e+00 J


## Step 35/796: t = 1.54e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 34 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.14e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.265e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.404e-04
  |ΔD|_∞ = 8.064e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.14e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.807e-05
  |ΔD|_∞ = 1.613e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.14e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.615e-06
  |ΔD|_∞ = 3.226e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.14e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.640e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.123e-06
  |ΔD|_∞ = 6.451e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.14e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.640e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.246e-07
  |ΔD|_∞ = 1.290e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.14e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.373e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.492e-08
  |ΔD|_∞ = 2.580e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4823e+00 J
  → Fracture energy : 1.4774e+00 J
  → Total energy    : 2.9596e+00 J


## Step 36/796: t = 1.58e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 35 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.15e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.260e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.321e-04
  |ΔD|_∞ = 7.078e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.15e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.384e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.641e-05
  |ΔD|_∞ = 1.416e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.15e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.384e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.283e-06
  |ΔD|_∞ = 2.831e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.15e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.598e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.057e-06
  |ΔD|_∞ = 5.663e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.15e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.598e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.113e-07
  |ΔD|_∞ = 1.133e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.15e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.906e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.226e-08
  |ΔD|_∞ = 2.265e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4858e+00 J
  → Fracture energy : 1.4775e+00 J
  → Total energy    : 2.9633e+00 J


## Step 37/796: t = 1.63e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 36 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.16e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.256e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.266e-04
  |ΔD|_∞ = 6.338e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.16e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.044e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.531e-05
  |ΔD|_∞ = 1.268e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.16e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.063e-06
  |ΔD|_∞ = 2.535e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.16e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.013e-06
  |ΔD|_∞ = 5.071e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.16e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.025e-07
  |ΔD|_∞ = 1.014e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.16e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.737e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.050e-08
  |ΔD|_∞ = 2.028e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4893e+00 J
  → Fracture energy : 1.4776e+00 J
  → Total energy    : 2.9669e+00 J


## Step 38/796: t = 1.68e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 37 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.17e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.253e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.229e-04
  |ΔD|_∞ = 5.786e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.17e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.177e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.458e-05
  |ΔD|_∞ = 1.157e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.17e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.177e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.916e-06
  |ΔD|_∞ = 2.314e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.17e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.739e-20

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.832e-07
  |ΔD|_∞ = 4.629e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.17e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.739e-20

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.966e-07
  |ΔD|_∞ = 9.257e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.17e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.692e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.933e-08
  |ΔD|_∞ = 1.851e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4928e+00 J
  → Fracture energy : 1.4777e+00 J
  → Total energy    : 2.9705e+00 J


## Step 39/796: t = 1.72e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 38 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.18e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.251e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.203e-04
  |ΔD|_∞ = 5.381e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.18e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.407e-05
  |ΔD|_∞ = 1.076e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.18e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.813e-06
  |ΔD|_∞ = 2.153e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.18e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.266e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.627e-07
  |ΔD|_∞ = 4.305e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.18e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.266e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.925e-07
  |ΔD|_∞ = 8.610e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.18e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.851e-08
  |ΔD|_∞ = 1.722e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4963e+00 J
  → Fracture energy : 1.4778e+00 J
  → Total energy    : 2.9742e+00 J


## Step 40/796: t = 1.77e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 39 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.19e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.248e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.186e-04
  |ΔD|_∞ = 5.090e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.19e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.373e-05
  |ΔD|_∞ = 1.018e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.19e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.745e-06
  |ΔD|_∞ = 2.036e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.19e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.490e-07
  |ΔD|_∞ = 4.072e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.19e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.559e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.898e-07
  |ΔD|_∞ = 8.144e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.19e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.559e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.796e-08
  |ΔD|_∞ = 1.629e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4999e+00 J
  → Fracture energy : 1.4779e+00 J
  → Total energy    : 2.9778e+00 J


## Step 41/796: t = 1.81e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 40 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.246e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.175e-04
  |ΔD|_∞ = 5.061e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.160e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.350e-05
  |ΔD|_∞ = 1.012e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.699e-06
  |ΔD|_∞ = 2.024e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.399e-07
  |ΔD|_∞ = 4.049e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.880e-07
  |ΔD|_∞ = 8.097e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.366e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.760e-08
  |ΔD|_∞ = 1.619e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5034e+00 J
  → Fracture energy : 1.4780e+00 J
  → Total energy    : 2.9814e+00 J


## Step 42/796: t = 1.86e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 41 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.21e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.244e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.167e-04
  |ΔD|_∞ = 5.045e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.21e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.846e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.334e-05
  |ΔD|_∞ = 1.009e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.21e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.846e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.668e-06
  |ΔD|_∞ = 2.018e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.21e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.337e-07
  |ΔD|_∞ = 4.036e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.21e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.867e-07
  |ΔD|_∞ = 8.072e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.21e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.735e-08
  |ΔD|_∞ = 1.614e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5070e+00 J
  → Fracture energy : 1.4781e+00 J
  → Total energy    : 2.9851e+00 J


## Step 43/796: t = 1.90e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 42 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.22e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.243e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.162e-04
  |ΔD|_∞ = 5.039e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.22e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.684e-20

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.324e-05
  |ΔD|_∞ = 1.008e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.22e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.030e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.648e-06
  |ΔD|_∞ = 2.016e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.22e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.296e-07
  |ΔD|_∞ = 4.031e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.22e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.385e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.859e-07
  |ΔD|_∞ = 8.063e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.22e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.385e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.719e-08
  |ΔD|_∞ = 1.613e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5105e+00 J
  → Fracture energy : 1.4782e+00 J
  → Total energy    : 2.9887e+00 J


## Step 44/796: t = 1.95e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 43 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.23e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.241e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.159e-04
  |ΔD|_∞ = 5.040e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.23e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.318e-05
  |ΔD|_∞ = 1.008e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.23e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.635e-06
  |ΔD|_∞ = 2.016e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.23e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.271e-07
  |ΔD|_∞ = 4.032e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.23e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.854e-07
  |ΔD|_∞ = 8.064e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.23e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.708e-08
  |ΔD|_∞ = 1.613e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5141e+00 J
  → Fracture energy : 1.4784e+00 J
  → Total energy    : 2.9924e+00 J


## Step 45/796: t = 1.99e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 44 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.24e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.239e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.157e-04
  |ΔD|_∞ = 5.046e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.24e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.709e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.314e-05
  |ΔD|_∞ = 1.009e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.24e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.709e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.628e-06
  |ΔD|_∞ = 2.018e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.24e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.709e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.255e-07
  |ΔD|_∞ = 4.037e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.24e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.851e-07
  |ΔD|_∞ = 8.073e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.24e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.709e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.702e-08
  |ΔD|_∞ = 1.615e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5176e+00 J
  → Fracture energy : 1.4785e+00 J
  → Total energy    : 2.9961e+00 J


## Step 46/796: t = 2.04e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 45 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.25e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.238e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.156e-04
  |ΔD|_∞ = 5.055e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.25e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.454e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.312e-05
  |ΔD|_∞ = 1.011e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.25e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.623e-06
  |ΔD|_∞ = 2.022e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.25e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.515e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.247e-07
  |ΔD|_∞ = 4.044e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.25e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.515e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.849e-07
  |ΔD|_∞ = 8.088e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.25e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.699e-08
  |ΔD|_∞ = 1.618e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5212e+00 J
  → Fracture energy : 1.4786e+00 J
  → Total energy    : 2.9997e+00 J


## Step 47/796: t = 2.08e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 46 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.26e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.236e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.155e-04
  |ΔD|_∞ = 5.067e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.26e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.311e-05
  |ΔD|_∞ = 1.013e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.26e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.622e-06
  |ΔD|_∞ = 2.027e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.26e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.243e-07
  |ΔD|_∞ = 4.054e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.26e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.849e-07
  |ΔD|_∞ = 8.107e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.26e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.769e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.697e-08
  |ΔD|_∞ = 1.621e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5247e+00 J
  → Fracture energy : 1.4787e+00 J
  → Total energy    : 3.0034e+00 J


## Step 48/796: t = 2.13e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 47 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.27e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.234e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.155e-04
  |ΔD|_∞ = 5.080e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.27e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.311e-05
  |ΔD|_∞ = 1.016e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.27e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.622e-06
  |ΔD|_∞ = 2.032e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.27e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.333e-20

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.244e-07
  |ΔD|_∞ = 4.064e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.27e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.333e-20

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.849e-07
  |ΔD|_∞ = 8.128e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.27e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.697e-08
  |ΔD|_∞ = 1.626e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5283e+00 J
  → Fracture energy : 1.4788e+00 J
  → Total energy    : 3.0071e+00 J


## Step 49/796: t = 2.17e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 48 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.28e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.233e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.156e-04
  |ΔD|_∞ = 5.095e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.28e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.311e-05
  |ΔD|_∞ = 1.019e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.28e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.623e-06
  |ΔD|_∞ = 2.038e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.28e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.246e-07
  |ΔD|_∞ = 4.076e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.28e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.849e-07
  |ΔD|_∞ = 8.152e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.28e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.684e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.698e-08
  |ΔD|_∞ = 1.630e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5319e+00 J
  → Fracture energy : 1.4789e+00 J
  → Total energy    : 3.0108e+00 J


## Step 50/796: t = 2.22e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 49 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.29e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.231e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.156e-04
  |ΔD|_∞ = 5.110e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.29e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.313e-05
  |ΔD|_∞ = 1.022e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.29e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.625e-06
  |ΔD|_∞ = 2.044e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.29e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.250e-07
  |ΔD|_∞ = 4.088e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.29e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.850e-07
  |ΔD|_∞ = 8.177e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.29e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.700e-08
  |ΔD|_∞ = 1.635e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5355e+00 J
  → Fracture energy : 1.4790e+00 J
  → Total energy    : 3.0144e+00 J


## Step 51/796: t = 2.26e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 50 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.3e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.230e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.156e-04
  |ΔD|_∞ = 5.124e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.3e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.313e-05
  |ΔD|_∞ = 1.025e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.3e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.625e-06
  |ΔD|_∞ = 2.050e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.3e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.251e-07
  |ΔD|_∞ = 4.099e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.3e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.850e-07
  |ΔD|_∞ = 8.198e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.3e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.700e-08
  |ΔD|_∞ = 1.640e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5390e+00 J
  → Fracture energy : 1.4791e+00 J
  → Total energy    : 3.0181e+00 J


## Step 52/796: t = 2.31e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 51 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.31e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.228e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.259e-04
  |ΔD|_∞ = 1.698e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.31e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.519e-05
  |ΔD|_∞ = 3.396e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.31e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.037e-06
  |ΔD|_∞ = 6.792e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.31e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.007e-06
  |ΔD|_∞ = 1.358e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.31e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.015e-07
  |ΔD|_∞ = 2.717e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.31e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.030e-08
  |ΔD|_∞ = 5.434e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5426e+00 J
  → Fracture energy : 1.4792e+00 J
  → Total energy    : 3.0218e+00 J


## Step 53/796: t = 2.35e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 52 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.32e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.228e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.229e-04
  |ΔD|_∞ = 7.635e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.32e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.458e-05
  |ΔD|_∞ = 1.527e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.32e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.917e-06
  |ΔD|_∞ = 3.054e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.32e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.833e-07
  |ΔD|_∞ = 6.108e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.32e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.967e-07
  |ΔD|_∞ = 1.222e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.32e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.073e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.933e-08
  |ΔD|_∞ = 2.443e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5462e+00 J
  → Fracture energy : 1.4793e+00 J
  → Total energy    : 3.0255e+00 J


## Step 54/796: t = 2.40e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 53 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.33e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.226e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.209e-04
  |ΔD|_∞ = 5.327e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.33e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.461e-20

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.418e-05
  |ΔD|_∞ = 1.065e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.33e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.835e-06
  |ΔD|_∞ = 2.131e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.33e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.671e-07
  |ΔD|_∞ = 4.262e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.33e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.934e-07
  |ΔD|_∞ = 8.523e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.33e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.868e-08
  |ΔD|_∞ = 1.705e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5498e+00 J
  → Fracture energy : 1.4794e+00 J
  → Total energy    : 3.0292e+00 J


## Step 55/796: t = 2.45e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 54 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.34e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.224e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.192e-04
  |ΔD|_∞ = 5.314e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.34e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.960e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.384e-05
  |ΔD|_∞ = 1.063e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.34e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.082e-20

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.769e-06
  |ΔD|_∞ = 2.126e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.34e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.537e-07
  |ΔD|_∞ = 4.252e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.34e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.082e-20

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.907e-07
  |ΔD|_∞ = 8.503e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.34e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.082e-20

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.815e-08
  |ΔD|_∞ = 1.701e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5534e+00 J
  → Fracture energy : 1.4795e+00 J
  → Total energy    : 3.0329e+00 J


## Step 56/796: t = 2.49e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 55 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.35e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.222e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.181e-04
  |ΔD|_∞ = 5.299e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.35e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.362e-05
  |ΔD|_∞ = 1.060e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.35e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.691e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.725e-06
  |ΔD|_∞ = 2.119e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.35e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.691e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.450e-07
  |ΔD|_∞ = 4.239e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.35e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.165e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.890e-07
  |ΔD|_∞ = 8.478e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.35e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.165e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.780e-08
  |ΔD|_∞ = 1.696e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5570e+00 J
  → Fracture energy : 1.4797e+00 J
  → Total energy    : 3.0366e+00 J


## Step 57/796: t = 2.54e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 56 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.36e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.221e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.174e-04
  |ΔD|_∞ = 5.288e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.36e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.347e-05
  |ΔD|_∞ = 1.058e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.36e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.695e-06
  |ΔD|_∞ = 2.115e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.36e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.389e-07
  |ΔD|_∞ = 4.231e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.36e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.878e-07
  |ΔD|_∞ = 8.461e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.36e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.756e-08
  |ΔD|_∞ = 1.692e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5606e+00 J
  → Fracture energy : 1.4798e+00 J
  → Total energy    : 3.0403e+00 J


## Step 58/796: t = 2.58e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 57 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.37e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.219e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.169e-04
  |ΔD|_∞ = 5.284e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.37e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.338e-05
  |ΔD|_∞ = 1.057e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.37e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.676e-06
  |ΔD|_∞ = 2.114e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.37e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.352e-07
  |ΔD|_∞ = 4.228e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.37e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.870e-07
  |ΔD|_∞ = 8.455e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.37e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.741e-08
  |ΔD|_∞ = 1.691e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5642e+00 J
  → Fracture energy : 1.4799e+00 J
  → Total energy    : 3.0441e+00 J


## Step 59/796: t = 2.63e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 58 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.38e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.217e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.166e-04
  |ΔD|_∞ = 5.287e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.38e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.332e-05
  |ΔD|_∞ = 1.057e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.38e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.011e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.664e-06
  |ΔD|_∞ = 2.115e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.38e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.011e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.328e-07
  |ΔD|_∞ = 4.229e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.38e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.866e-07
  |ΔD|_∞ = 8.459e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.38e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.731e-08
  |ΔD|_∞ = 1.692e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5678e+00 J
  → Fracture energy : 1.4800e+00 J
  → Total energy    : 3.0478e+00 J


## Step 60/796: t = 2.67e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 59 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.39e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.216e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.166e-04
  |ΔD|_∞ = 5.293e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.39e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.331e-05
  |ΔD|_∞ = 1.059e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.39e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.662e-06
  |ΔD|_∞ = 2.117e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.39e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.324e-07
  |ΔD|_∞ = 4.235e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.39e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.865e-07
  |ΔD|_∞ = 8.469e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.39e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.730e-08
  |ΔD|_∞ = 1.694e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5714e+00 J
  → Fracture energy : 1.4801e+00 J
  → Total energy    : 3.0515e+00 J


## Step 61/796: t = 2.72e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 60 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.214e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.163e-04
  |ΔD|_∞ = 5.303e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.327e-05
  |ΔD|_∞ = 1.061e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.653e-06
  |ΔD|_∞ = 2.121e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.306e-07
  |ΔD|_∞ = 4.242e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.861e-07
  |ΔD|_∞ = 8.484e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.688e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.722e-08
  |ΔD|_∞ = 1.697e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5751e+00 J
  → Fracture energy : 1.4802e+00 J
  → Total energy    : 3.0552e+00 J


## Step 62/796: t = 2.76e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 61 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.41e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.213e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.163e-04
  |ΔD|_∞ = 5.315e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.41e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.487e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.326e-05
  |ΔD|_∞ = 1.063e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.41e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.651e-06
  |ΔD|_∞ = 2.126e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.41e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.709e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.303e-07
  |ΔD|_∞ = 4.252e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.41e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.709e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.861e-07
  |ΔD|_∞ = 8.504e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.41e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.721e-08
  |ΔD|_∞ = 1.701e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5787e+00 J
  → Fracture energy : 1.4803e+00 J
  → Total energy    : 3.0590e+00 J


## Step 63/796: t = 2.81e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 62 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.42e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.211e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.163e-04
  |ΔD|_∞ = 5.328e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.42e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.326e-05
  |ΔD|_∞ = 1.066e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.42e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.428e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.651e-06
  |ΔD|_∞ = 2.131e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.42e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.428e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.302e-07
  |ΔD|_∞ = 4.263e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.42e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.860e-07
  |ΔD|_∞ = 8.525e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.42e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.721e-08
  |ΔD|_∞ = 1.705e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5823e+00 J
  → Fracture energy : 1.4804e+00 J
  → Total energy    : 3.0627e+00 J


## Step 64/796: t = 2.85e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 63 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.43e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.210e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.199e-04
  |ΔD|_∞ = 1.685e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.43e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.643e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.397e-05
  |ΔD|_∞ = 3.370e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.43e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.794e-06
  |ΔD|_∞ = 6.739e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.43e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.588e-07
  |ΔD|_∞ = 1.348e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.43e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.918e-07
  |ΔD|_∞ = 2.696e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.43e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.643e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.835e-08
  |ΔD|_∞ = 5.391e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5859e+00 J
  → Fracture energy : 1.4805e+00 J
  → Total energy    : 3.0665e+00 J


## Step 65/796: t = 2.90e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 64 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.44e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.210e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.170e-04
  |ΔD|_∞ = 5.374e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.44e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.340e-05
  |ΔD|_∞ = 1.075e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.44e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.680e-06
  |ΔD|_∞ = 2.150e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.44e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.361e-07
  |ΔD|_∞ = 4.300e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.44e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.872e-07
  |ΔD|_∞ = 8.599e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.44e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.744e-08
  |ΔD|_∞ = 1.720e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5896e+00 J
  → Fracture energy : 1.4806e+00 J
  → Total energy    : 3.0702e+00 J


## Step 66/796: t = 2.94e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 65 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.45e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.207e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.168e-04
  |ΔD|_∞ = 5.386e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.45e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.335e-05
  |ΔD|_∞ = 1.077e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.45e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.671e-06
  |ΔD|_∞ = 2.155e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.45e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.341e-07
  |ΔD|_∞ = 4.309e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.45e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.868e-07
  |ΔD|_∞ = 8.618e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.45e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.736e-08
  |ΔD|_∞ = 1.724e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5932e+00 J
  → Fracture energy : 1.4807e+00 J
  → Total energy    : 3.0740e+00 J


## Step 67/796: t = 2.99e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 66 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.46e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.205e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.195e-04
  |ΔD|_∞ = 8.539e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.46e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.391e-05
  |ΔD|_∞ = 1.708e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.46e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.782e-06
  |ΔD|_∞ = 3.416e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.46e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.564e-07
  |ΔD|_∞ = 6.832e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.46e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.913e-07
  |ΔD|_∞ = 1.366e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.46e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.825e-08
  |ΔD|_∞ = 2.733e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5969e+00 J
  → Fracture energy : 1.4808e+00 J
  → Total energy    : 3.0777e+00 J


## Step 68/796: t = 3.03e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 67 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.47e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.204e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.192e-04
  |ΔD|_∞ = 5.480e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.47e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.384e-05
  |ΔD|_∞ = 1.096e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.47e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.768e-06
  |ΔD|_∞ = 2.192e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.47e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.536e-07
  |ΔD|_∞ = 4.384e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.47e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.907e-07
  |ΔD|_∞ = 8.768e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.47e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.814e-08
  |ΔD|_∞ = 1.754e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6005e+00 J
  → Fracture energy : 1.4810e+00 J
  → Total energy    : 3.0815e+00 J


## Step 69/796: t = 3.08e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 68 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.48e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.203e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.185e-04
  |ΔD|_∞ = 5.496e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.48e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.370e-05
  |ΔD|_∞ = 1.099e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.48e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.068e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.740e-06
  |ΔD|_∞ = 2.198e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.48e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.068e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.480e-07
  |ΔD|_∞ = 4.396e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.48e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.896e-07
  |ΔD|_∞ = 8.793e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.48e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.792e-08
  |ΔD|_∞ = 1.759e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6042e+00 J
  → Fracture energy : 1.4811e+00 J
  → Total energy    : 3.0852e+00 J


## Step 70/796: t = 3.12e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 69 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.49e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.201e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.179e-04
  |ΔD|_∞ = 5.498e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.49e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.359e-05
  |ΔD|_∞ = 1.100e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.49e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.718e-06
  |ΔD|_∞ = 2.199e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.49e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.436e-07
  |ΔD|_∞ = 4.398e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.49e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.887e-07
  |ΔD|_∞ = 8.797e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.49e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.774e-08
  |ΔD|_∞ = 1.759e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6078e+00 J
  → Fracture energy : 1.4812e+00 J
  → Total energy    : 3.0890e+00 J


## Step 71/796: t = 3.17e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 70 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.5e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.199e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.228e-04
  |ΔD|_∞ = 2.245e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.5e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.456e-05
  |ΔD|_∞ = 4.490e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.5e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.913e-06
  |ΔD|_∞ = 8.980e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.5e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.826e-07
  |ΔD|_∞ = 1.796e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.5e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.965e-07
  |ΔD|_∞ = 3.592e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.5e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.930e-08
  |ΔD|_∞ = 7.184e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6115e+00 J
  → Fracture energy : 1.4813e+00 J
  → Total energy    : 3.0928e+00 J


## Step 72/796: t = 3.22e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 71 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.51e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.204e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.184e-04
  |ΔD|_∞ = 5.534e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.51e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.368e-05
  |ΔD|_∞ = 1.107e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.51e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.735e-06
  |ΔD|_∞ = 2.213e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.51e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.843e-20

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.470e-07
  |ΔD|_∞ = 4.427e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.51e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.843e-20

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.894e-07
  |ΔD|_∞ = 8.854e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.51e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.788e-08
  |ΔD|_∞ = 1.771e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6151e+00 J
  → Fracture energy : 1.4814e+00 J
  → Total energy    : 3.0965e+00 J


## Step 73/796: t = 3.26e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 72 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.52e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.198e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.186e-04
  |ΔD|_∞ = 5.535e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.52e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.371e-05
  |ΔD|_∞ = 1.107e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.52e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.742e-06
  |ΔD|_∞ = 2.214e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.52e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.485e-07
  |ΔD|_∞ = 4.428e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.52e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.897e-07
  |ΔD|_∞ = 8.856e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.52e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.794e-08
  |ΔD|_∞ = 1.771e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6188e+00 J
  → Fracture energy : 1.4815e+00 J
  → Total energy    : 3.1003e+00 J


## Step 74/796: t = 3.31e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 73 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.53e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.195e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.183e-04
  |ΔD|_∞ = 5.563e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.53e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.367e-05
  |ΔD|_∞ = 1.113e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.53e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.733e-06
  |ΔD|_∞ = 2.225e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.53e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.466e-07
  |ΔD|_∞ = 4.451e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.53e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.893e-07
  |ΔD|_∞ = 8.901e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.53e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.787e-08
  |ΔD|_∞ = 1.780e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6225e+00 J
  → Fracture energy : 1.4816e+00 J
  → Total energy    : 3.1041e+00 J


## Step 75/796: t = 3.35e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 74 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.54e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.194e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.182e-04
  |ΔD|_∞ = 5.571e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.54e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.365e-05
  |ΔD|_∞ = 1.114e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.54e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.729e-06
  |ΔD|_∞ = 2.228e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.54e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.986e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.459e-07
  |ΔD|_∞ = 4.457e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.54e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.986e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.892e-07
  |ΔD|_∞ = 8.914e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.54e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.784e-08
  |ΔD|_∞ = 1.783e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6262e+00 J
  → Fracture energy : 1.4817e+00 J
  → Total energy    : 3.1079e+00 J


## Step 76/796: t = 3.40e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 75 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.55e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.192e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.181e-04
  |ΔD|_∞ = 5.586e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.55e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.363e-05
  |ΔD|_∞ = 1.117e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.55e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.726e-06
  |ΔD|_∞ = 2.235e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.55e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.452e-07
  |ΔD|_∞ = 4.469e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.55e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.890e-07
  |ΔD|_∞ = 8.938e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.55e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.781e-08
  |ΔD|_∞ = 1.788e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6298e+00 J
  → Fracture energy : 1.4818e+00 J
  → Total energy    : 3.1117e+00 J


## Step 77/796: t = 3.44e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 76 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.56e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.191e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.181e-04
  |ΔD|_∞ = 5.599e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.56e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.418e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.361e-05
  |ΔD|_∞ = 1.120e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.56e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.723e-06
  |ΔD|_∞ = 2.240e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.56e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.446e-07
  |ΔD|_∞ = 4.480e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.56e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.889e-07
  |ΔD|_∞ = 8.959e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.56e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.778e-08
  |ΔD|_∞ = 1.792e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6335e+00 J
  → Fracture energy : 1.4820e+00 J
  → Total energy    : 3.1155e+00 J


## Step 78/796: t = 3.49e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 77 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.57e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.189e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.180e-04
  |ΔD|_∞ = 5.613e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.57e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.361e-05
  |ΔD|_∞ = 1.123e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.57e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.721e-06
  |ΔD|_∞ = 2.245e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.57e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.443e-07
  |ΔD|_∞ = 4.491e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.57e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.050e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.889e-07
  |ΔD|_∞ = 8.981e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.57e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.954e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.777e-08
  |ΔD|_∞ = 1.796e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6372e+00 J
  → Fracture energy : 1.4821e+00 J
  → Total energy    : 3.1193e+00 J


## Step 79/796: t = 3.53e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 78 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.58e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.188e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.180e-04
  |ΔD|_∞ = 5.629e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.58e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.361e-05
  |ΔD|_∞ = 1.126e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.58e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.722e-06
  |ΔD|_∞ = 2.251e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.58e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.444e-07
  |ΔD|_∞ = 4.503e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.58e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.889e-07
  |ΔD|_∞ = 9.006e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.58e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.777e-08
  |ΔD|_∞ = 1.801e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6409e+00 J
  → Fracture energy : 1.4822e+00 J
  → Total energy    : 3.1231e+00 J


## Step 80/796: t = 3.58e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 79 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.59e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.186e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.181e-04
  |ΔD|_∞ = 5.645e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.59e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.754e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.362e-05
  |ΔD|_∞ = 1.129e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.59e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.724e-06
  |ΔD|_∞ = 2.258e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.59e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.447e-07
  |ΔD|_∞ = 4.516e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.59e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.754e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.889e-07
  |ΔD|_∞ = 9.032e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.59e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.754e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.779e-08
  |ΔD|_∞ = 1.806e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6446e+00 J
  → Fracture energy : 1.4823e+00 J
  → Total energy    : 3.1269e+00 J


## Step 81/796: t = 3.62e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 80 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.185e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.183e-04
  |ΔD|_∞ = 5.662e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.367e-05
  |ΔD|_∞ = 1.132e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.733e-06
  |ΔD|_∞ = 2.265e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.466e-07
  |ΔD|_∞ = 4.530e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.893e-07
  |ΔD|_∞ = 9.060e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.562e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.787e-08
  |ΔD|_∞ = 1.812e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6483e+00 J
  → Fracture energy : 1.4824e+00 J
  → Total energy    : 3.1307e+00 J


## Step 82/796: t = 3.67e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 81 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.61e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.183e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.182e-04
  |ΔD|_∞ = 5.680e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.61e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.365e-05
  |ΔD|_∞ = 1.136e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.61e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.729e-06
  |ΔD|_∞ = 2.272e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.61e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.459e-07
  |ΔD|_∞ = 4.544e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.61e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.645e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.892e-07
  |ΔD|_∞ = 9.089e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.61e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.836e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.783e-08
  |ΔD|_∞ = 1.818e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6520e+00 J
  → Fracture energy : 1.4825e+00 J
  → Total energy    : 3.1345e+00 J


## Step 83/796: t = 3.71e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 82 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.62e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.182e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.183e-04
  |ΔD|_∞ = 5.699e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.62e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.367e-05
  |ΔD|_∞ = 1.140e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.62e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.463e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.733e-06
  |ΔD|_∞ = 2.280e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.62e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.463e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.466e-07
  |ΔD|_∞ = 4.559e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.62e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.893e-07
  |ΔD|_∞ = 9.119e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.62e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.787e-08
  |ΔD|_∞ = 1.824e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6557e+00 J
  → Fracture energy : 1.4826e+00 J
  → Total energy    : 3.1383e+00 J


## Step 84/796: t = 3.76e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 83 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.63e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.181e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.184e-04
  |ΔD|_∞ = 5.718e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.63e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.019e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.369e-05
  |ΔD|_∞ = 1.144e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.63e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.019e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.737e-06
  |ΔD|_∞ = 2.287e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.63e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.475e-07
  |ΔD|_∞ = 4.575e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.63e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.895e-07
  |ΔD|_∞ = 9.149e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.63e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.790e-08
  |ΔD|_∞ = 1.830e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6594e+00 J
  → Fracture energy : 1.4827e+00 J
  → Total energy    : 3.1422e+00 J


## Step 85/796: t = 3.80e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 84 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.64e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.179e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.185e-04
  |ΔD|_∞ = 5.738e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.64e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.040e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.371e-05
  |ΔD|_∞ = 1.148e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.64e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.742e-06
  |ΔD|_∞ = 2.295e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.64e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.483e-07
  |ΔD|_∞ = 4.590e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.64e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.897e-07
  |ΔD|_∞ = 9.180e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.64e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.793e-08
  |ΔD|_∞ = 1.836e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6632e+00 J
  → Fracture energy : 1.4828e+00 J
  → Total energy    : 3.1460e+00 J


## Step 86/796: t = 3.85e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 85 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.65e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.178e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.187e-04
  |ΔD|_∞ = 5.757e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.65e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.644e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.373e-05
  |ΔD|_∞ = 1.151e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.65e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.611e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.746e-06
  |ΔD|_∞ = 2.303e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.65e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.492e-07
  |ΔD|_∞ = 4.606e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.65e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.898e-07
  |ΔD|_∞ = 9.212e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.65e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.367e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.797e-08
  |ΔD|_∞ = 1.842e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6669e+00 J
  → Fracture energy : 1.4829e+00 J
  → Total energy    : 3.1498e+00 J


## Step 87/796: t = 3.89e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 86 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.66e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.176e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.188e-04
  |ΔD|_∞ = 5.777e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.66e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.375e-05
  |ΔD|_∞ = 1.155e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.66e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.751e-06
  |ΔD|_∞ = 2.311e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.66e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.502e-07
  |ΔD|_∞ = 4.622e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.66e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.900e-07
  |ΔD|_∞ = 9.244e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.66e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.801e-08
  |ΔD|_∞ = 1.849e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6706e+00 J
  → Fracture energy : 1.4831e+00 J
  → Total energy    : 3.1537e+00 J


## Step 88/796: t = 3.94e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 87 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.67e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.175e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.189e-04
  |ΔD|_∞ = 5.798e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.67e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.378e-05
  |ΔD|_∞ = 1.160e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.67e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.756e-06
  |ΔD|_∞ = 2.319e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.67e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.453e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.512e-07
  |ΔD|_∞ = 4.638e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.67e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.453e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.902e-07
  |ΔD|_∞ = 9.276e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.67e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.293e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.805e-08
  |ΔD|_∞ = 1.855e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6743e+00 J
  → Fracture energy : 1.4832e+00 J
  → Total energy    : 3.1575e+00 J


## Step 89/796: t = 3.98e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 88 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.68e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.174e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.190e-04
  |ΔD|_∞ = 5.818e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.68e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.381e-05
  |ΔD|_∞ = 1.164e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.68e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.761e-06
  |ΔD|_∞ = 2.327e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.68e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.522e-07
  |ΔD|_∞ = 4.654e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.68e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.904e-07
  |ΔD|_∞ = 9.309e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.68e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.809e-08
  |ΔD|_∞ = 1.862e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6781e+00 J
  → Fracture energy : 1.4833e+00 J
  → Total energy    : 3.1614e+00 J


## Step 90/796: t = 4.03e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 89 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.69e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.172e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.192e-04
  |ΔD|_∞ = 5.838e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.69e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.611e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.384e-05
  |ΔD|_∞ = 1.168e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.69e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.004e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.769e-06
  |ΔD|_∞ = 2.335e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.69e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.538e-07
  |ΔD|_∞ = 4.671e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.69e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.908e-07
  |ΔD|_∞ = 9.341e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.69e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.815e-08
  |ΔD|_∞ = 1.868e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6818e+00 J
  → Fracture energy : 1.4834e+00 J
  → Total energy    : 3.1652e+00 J


## Step 91/796: t = 4.08e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 90 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.7e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.171e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.193e-04
  |ΔD|_∞ = 5.859e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.7e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.386e-05
  |ΔD|_∞ = 1.172e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.7e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.772e-06
  |ΔD|_∞ = 2.344e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.7e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.544e-07
  |ΔD|_∞ = 4.687e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.7e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.909e-07
  |ΔD|_∞ = 9.375e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.7e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.817e-08
  |ΔD|_∞ = 1.875e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6856e+00 J
  → Fracture energy : 1.4835e+00 J
  → Total energy    : 3.1691e+00 J


## Step 92/796: t = 4.12e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 91 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.71e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.169e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.194e-04
  |ΔD|_∞ = 5.880e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.71e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.389e-05
  |ΔD|_∞ = 1.176e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.71e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.777e-06
  |ΔD|_∞ = 2.352e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.71e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.554e-07
  |ΔD|_∞ = 4.704e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.71e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.911e-07
  |ΔD|_∞ = 9.408e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.71e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.822e-08
  |ΔD|_∞ = 1.882e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6893e+00 J
  → Fracture energy : 1.4836e+00 J
  → Total energy    : 3.1729e+00 J


## Step 93/796: t = 4.17e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 92 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.72e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.168e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.196e-04
  |ΔD|_∞ = 5.901e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.72e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.392e-05
  |ΔD|_∞ = 1.180e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.72e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.783e-06
  |ΔD|_∞ = 2.360e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.72e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.566e-07
  |ΔD|_∞ = 4.721e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.72e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.376e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.913e-07
  |ΔD|_∞ = 9.442e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.72e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.376e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.827e-08
  |ΔD|_∞ = 1.888e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6931e+00 J
  → Fracture energy : 1.4837e+00 J
  → Total energy    : 3.1768e+00 J


## Step 94/796: t = 4.21e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 93 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.73e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.167e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.197e-04
  |ΔD|_∞ = 5.922e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.73e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.394e-05
  |ΔD|_∞ = 1.184e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.73e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.788e-06
  |ΔD|_∞ = 2.369e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.73e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.505e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.577e-07
  |ΔD|_∞ = 4.738e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.73e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.505e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.915e-07
  |ΔD|_∞ = 9.476e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.73e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.831e-08
  |ΔD|_∞ = 1.895e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6968e+00 J
  → Fracture energy : 1.4838e+00 J
  → Total energy    : 3.1807e+00 J


## Step 95/796: t = 4.26e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 94 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.74e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.165e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.199e-04
  |ΔD|_∞ = 5.944e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.74e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.397e-05
  |ΔD|_∞ = 1.189e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.74e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.794e-06
  |ΔD|_∞ = 2.378e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.74e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.589e-07
  |ΔD|_∞ = 4.755e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.74e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.918e-07
  |ΔD|_∞ = 9.510e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.74e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.835e-08
  |ΔD|_∞ = 1.902e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7006e+00 J
  → Fracture energy : 1.4840e+00 J
  → Total energy    : 3.1845e+00 J


## Step 96/796: t = 4.30e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 95 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.75e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.164e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.200e-04
  |ΔD|_∞ = 5.965e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.75e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.421e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.400e-05
  |ΔD|_∞ = 1.193e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.75e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.800e-06
  |ΔD|_∞ = 2.386e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.75e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.488e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.601e-07
  |ΔD|_∞ = 4.772e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.75e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.488e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.920e-07
  |ΔD|_∞ = 9.545e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.75e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.840e-08
  |ΔD|_∞ = 1.909e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7043e+00 J
  → Fracture energy : 1.4841e+00 J
  → Total energy    : 3.1884e+00 J


## Step 97/796: t = 4.35e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 96 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.76e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.162e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.202e-04
  |ΔD|_∞ = 5.987e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.76e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.403e-05
  |ΔD|_∞ = 1.197e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.76e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.806e-06
  |ΔD|_∞ = 2.395e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.76e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.613e-07
  |ΔD|_∞ = 4.790e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.76e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.923e-07
  |ΔD|_∞ = 9.580e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.76e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.845e-08
  |ΔD|_∞ = 1.916e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7081e+00 J
  → Fracture energy : 1.4842e+00 J
  → Total energy    : 3.1923e+00 J


## Step 98/796: t = 4.39e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 97 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.77e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.161e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.203e-04
  |ΔD|_∞ = 6.009e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.77e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.406e-05
  |ΔD|_∞ = 1.202e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.77e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.813e-06
  |ΔD|_∞ = 2.404e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.77e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.626e-07
  |ΔD|_∞ = 4.808e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.77e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.925e-07
  |ΔD|_∞ = 9.615e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.77e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.850e-08
  |ΔD|_∞ = 1.923e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7119e+00 J
  → Fracture energy : 1.4843e+00 J
  → Total energy    : 3.1962e+00 J


## Step 99/796: t = 4.44e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 98 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.78e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.160e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.206e-04
  |ΔD|_∞ = 6.032e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.78e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.411e-05
  |ΔD|_∞ = 1.206e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.78e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.823e-06
  |ΔD|_∞ = 2.413e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.78e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.646e-07
  |ΔD|_∞ = 4.825e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.78e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.199e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.929e-07
  |ΔD|_∞ = 9.651e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.78e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.199e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.858e-08
  |ΔD|_∞ = 1.930e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7157e+00 J
  → Fracture energy : 1.4844e+00 J
  → Total energy    : 3.2001e+00 J


## Step 100/796: t = 4.48e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 99 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.79e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.158e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.206e-04
  |ΔD|_∞ = 6.054e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.79e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.413e-05
  |ΔD|_∞ = 1.211e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.79e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.826e-06
  |ΔD|_∞ = 2.422e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.79e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.651e-07
  |ΔD|_∞ = 4.843e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.79e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.930e-07
  |ΔD|_∞ = 9.686e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.79e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.861e-08
  |ΔD|_∞ = 1.937e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7194e+00 J
  → Fracture energy : 1.4845e+00 J
  → Total energy    : 3.2040e+00 J


## Step 101/796: t = 4.53e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 100 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.157e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.208e-04
  |ΔD|_∞ = 6.077e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.417e-05
  |ΔD|_∞ = 1.215e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.834e-06
  |ΔD|_∞ = 2.431e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.880e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.668e-07
  |ΔD|_∞ = 4.861e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.919e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.934e-07
  |ΔD|_∞ = 9.723e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.723e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.867e-08
  |ΔD|_∞ = 1.945e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7232e+00 J
  → Fracture energy : 1.4846e+00 J
  → Total energy    : 3.2079e+00 J


## Step 102/796: t = 4.57e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 101 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.81e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.155e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.210e-04
  |ΔD|_∞ = 6.100e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.81e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.677e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.420e-05
  |ΔD|_∞ = 1.220e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.81e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.839e-06
  |ΔD|_∞ = 2.440e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.81e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.679e-07
  |ΔD|_∞ = 4.880e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.81e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.936e-07
  |ΔD|_∞ = 9.759e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.81e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.871e-08
  |ΔD|_∞ = 1.952e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7270e+00 J
  → Fracture energy : 1.4847e+00 J
  → Total energy    : 3.2118e+00 J


## Step 103/796: t = 4.62e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 102 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.82e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.154e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.212e-04
  |ΔD|_∞ = 6.123e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.82e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.238e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.423e-05
  |ΔD|_∞ = 1.225e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.82e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.846e-06
  |ΔD|_∞ = 2.449e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.82e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.693e-07
  |ΔD|_∞ = 4.898e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.82e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.939e-07
  |ΔD|_∞ = 9.796e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.82e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.877e-08
  |ΔD|_∞ = 1.959e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7308e+00 J
  → Fracture energy : 1.4849e+00 J
  → Total energy    : 3.2157e+00 J


## Step 104/796: t = 4.66e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 103 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.83e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.153e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.213e-04
  |ΔD|_∞ = 6.146e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.83e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.427e-05
  |ΔD|_∞ = 1.229e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.83e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.854e-06
  |ΔD|_∞ = 2.458e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.83e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.707e-07
  |ΔD|_∞ = 4.917e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.83e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.941e-07
  |ΔD|_∞ = 9.834e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.83e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.883e-08
  |ΔD|_∞ = 1.967e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7346e+00 J
  → Fracture energy : 1.4850e+00 J
  → Total energy    : 3.2196e+00 J


## Step 105/796: t = 4.71e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 104 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.84e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.151e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.215e-04
  |ΔD|_∞ = 6.170e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.84e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.430e-05
  |ΔD|_∞ = 1.234e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.84e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.860e-06
  |ΔD|_∞ = 2.468e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.84e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.721e-07
  |ΔD|_∞ = 4.936e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.84e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.944e-07
  |ΔD|_∞ = 9.871e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.84e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.888e-08
  |ΔD|_∞ = 1.974e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7384e+00 J
  → Fracture energy : 1.4851e+00 J
  → Total energy    : 3.2235e+00 J


## Step 106/796: t = 4.75e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 105 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.85e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.150e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.219e-04
  |ΔD|_∞ = 6.193e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.85e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.438e-05
  |ΔD|_∞ = 1.239e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.85e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.875e-06
  |ΔD|_∞ = 2.477e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.85e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.750e-07
  |ΔD|_∞ = 4.955e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.85e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.950e-07
  |ΔD|_∞ = 9.909e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.85e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.900e-08
  |ΔD|_∞ = 1.982e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7422e+00 J
  → Fracture energy : 1.4852e+00 J
  → Total energy    : 3.2274e+00 J


## Step 107/796: t = 4.80e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 106 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.86e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.149e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.219e-04
  |ΔD|_∞ = 6.217e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.86e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.438e-05
  |ΔD|_∞ = 1.243e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.86e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.875e-06
  |ΔD|_∞ = 2.487e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.86e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.750e-07
  |ΔD|_∞ = 4.974e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.86e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.950e-07
  |ΔD|_∞ = 9.947e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.86e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.900e-08
  |ΔD|_∞ = 1.989e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7460e+00 J
  → Fracture energy : 1.4853e+00 J
  → Total energy    : 3.2313e+00 J


## Step 108/796: t = 4.85e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 107 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.87e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.147e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.221e-04
  |ΔD|_∞ = 6.241e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.87e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.441e-05
  |ΔD|_∞ = 1.248e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.87e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.883e-06
  |ΔD|_∞ = 2.497e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.87e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.766e-07
  |ΔD|_∞ = 4.993e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.87e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.953e-07
  |ΔD|_∞ = 9.986e-07

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.87e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.906e-08
  |ΔD|_∞ = 1.997e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7498e+00 J
  → Fracture energy : 1.4854e+00 J
  → Total energy    : 3.2352e+00 J


## Step 109/796: t = 4.89e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 108 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.88e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.146e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.223e-04
  |ΔD|_∞ = 6.266e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.88e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.445e-05
  |ΔD|_∞ = 1.253e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.88e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.891e-06
  |ΔD|_∞ = 2.506e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.88e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.782e-07
  |ΔD|_∞ = 5.013e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.88e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.919e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.956e-07
  |ΔD|_∞ = 1.003e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.88e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.563e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.913e-08
  |ΔD|_∞ = 2.005e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7536e+00 J
  → Fracture energy : 1.4855e+00 J
  → Total energy    : 3.2392e+00 J


## Step 110/796: t = 4.94e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 109 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.89e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.145e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.225e-04
  |ΔD|_∞ = 6.290e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.89e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.317e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.449e-05
  |ΔD|_∞ = 1.258e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.89e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.899e-06
  |ΔD|_∞ = 2.516e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.89e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.503e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.797e-07
  |ΔD|_∞ = 5.032e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.89e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.503e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.959e-07
  |ΔD|_∞ = 1.006e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.89e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.919e-08
  |ΔD|_∞ = 2.013e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7574e+00 J
  → Fracture energy : 1.4857e+00 J
  → Total energy    : 3.2431e+00 J


## Step 111/796: t = 4.98e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 110 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.9e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.143e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.227e-04
  |ΔD|_∞ = 6.315e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.9e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.750e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.453e-05
  |ΔD|_∞ = 1.263e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.9e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.750e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.906e-06
  |ΔD|_∞ = 2.526e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.9e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.813e-07
  |ΔD|_∞ = 5.052e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.9e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.963e-07
  |ΔD|_∞ = 1.010e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.9e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.925e-08
  |ΔD|_∞ = 2.021e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7613e+00 J
  → Fracture energy : 1.4858e+00 J
  → Total energy    : 3.2470e+00 J


## Step 112/796: t = 5.03e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 111 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.91e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.142e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.229e-04
  |ΔD|_∞ = 6.340e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.91e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.289e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.457e-05
  |ΔD|_∞ = 1.268e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.91e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.289e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.914e-06
  |ΔD|_∞ = 2.536e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.91e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.289e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.829e-07
  |ΔD|_∞ = 5.072e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.91e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.966e-07
  |ΔD|_∞ = 1.014e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.91e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.932e-08
  |ΔD|_∞ = 2.029e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7651e+00 J
  → Fracture energy : 1.4859e+00 J
  → Total energy    : 3.2510e+00 J


## Step 113/796: t = 5.07e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 112 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.92e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.141e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.231e-04
  |ΔD|_∞ = 6.366e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.92e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.461e-05
  |ΔD|_∞ = 1.273e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.92e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.923e-06
  |ΔD|_∞ = 2.546e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.92e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.845e-07
  |ΔD|_∞ = 5.092e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.92e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.969e-07
  |ΔD|_∞ = 1.018e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.92e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.938e-08
  |ΔD|_∞ = 2.037e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7689e+00 J
  → Fracture energy : 1.4860e+00 J
  → Total energy    : 3.2549e+00 J


## Step 114/796: t = 5.12e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 113 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.93e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.139e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.233e-04
  |ΔD|_∞ = 6.391e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.93e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.466e-05
  |ΔD|_∞ = 1.278e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.93e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.931e-06
  |ΔD|_∞ = 2.556e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.93e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.862e-07
  |ΔD|_∞ = 5.113e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.93e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.972e-07
  |ΔD|_∞ = 1.023e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.93e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.945e-08
  |ΔD|_∞ = 2.045e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7728e+00 J
  → Fracture energy : 1.4861e+00 J
  → Total energy    : 3.2589e+00 J


## Step 115/796: t = 5.16e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 114 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.94e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.138e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.236e-04
  |ΔD|_∞ = 6.417e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.94e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.471e-05
  |ΔD|_∞ = 1.283e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.94e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.128e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.943e-06
  |ΔD|_∞ = 2.567e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.94e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.128e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.885e-07
  |ΔD|_∞ = 5.134e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.94e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.977e-07
  |ΔD|_∞ = 1.027e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.94e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.954e-08
  |ΔD|_∞ = 2.053e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7766e+00 J
  → Fracture energy : 1.4862e+00 J
  → Total energy    : 3.2628e+00 J


## Step 116/796: t = 5.21e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 115 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.95e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.136e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.237e-04
  |ΔD|_∞ = 6.443e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.95e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.815e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.474e-05
  |ΔD|_∞ = 1.289e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.95e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.842e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.948e-06
  |ΔD|_∞ = 2.577e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.95e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.842e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.897e-07
  |ΔD|_∞ = 5.154e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.95e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.979e-07
  |ΔD|_∞ = 1.031e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.95e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.959e-08
  |ΔD|_∞ = 2.062e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7804e+00 J
  → Fracture energy : 1.4864e+00 J
  → Total energy    : 3.2668e+00 J


## Step 117/796: t = 5.25e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 116 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.96e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.135e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.239e-04
  |ΔD|_∞ = 6.469e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.96e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.479e-05
  |ΔD|_∞ = 1.294e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.96e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.957e-06
  |ΔD|_∞ = 2.588e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.96e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.539e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.915e-07
  |ΔD|_∞ = 5.175e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.96e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.539e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.983e-07
  |ΔD|_∞ = 1.035e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.96e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.966e-08
  |ΔD|_∞ = 2.070e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7843e+00 J
  → Fracture energy : 1.4865e+00 J
  → Total energy    : 3.2708e+00 J


## Step 118/796: t = 5.30e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 117 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.97e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.134e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.242e-04
  |ΔD|_∞ = 6.496e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.97e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.559e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.483e-05
  |ΔD|_∞ = 1.299e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.97e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.673e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.966e-06
  |ΔD|_∞ = 2.598e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.97e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.673e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.933e-07
  |ΔD|_∞ = 5.197e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.97e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.987e-07
  |ΔD|_∞ = 1.039e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.97e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.973e-08
  |ΔD|_∞ = 2.079e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7881e+00 J
  → Fracture energy : 1.4866e+00 J
  → Total energy    : 3.2747e+00 J


## Step 119/796: t = 5.34e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 118 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.98e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.133e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.244e-04
  |ΔD|_∞ = 6.523e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.98e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.488e-05
  |ΔD|_∞ = 1.305e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.98e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.976e-06
  |ΔD|_∞ = 2.609e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.98e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.951e-07
  |ΔD|_∞ = 5.218e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.98e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.990e-07
  |ΔD|_∞ = 1.044e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.98e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.980e-08
  |ΔD|_∞ = 2.087e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7920e+00 J
  → Fracture energy : 1.4867e+00 J
  → Total energy    : 3.2787e+00 J


## Step 120/796: t = 5.39e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 119 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.99e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.131e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.246e-04
  |ΔD|_∞ = 6.550e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.99e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.464e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.492e-05
  |ΔD|_∞ = 1.310e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.99e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.464e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.985e-06
  |ΔD|_∞ = 2.620e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.99e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.970e-07
  |ΔD|_∞ = 5.240e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.99e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.994e-07
  |ΔD|_∞ = 1.048e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.99e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.988e-08
  |ΔD|_∞ = 2.096e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7958e+00 J
  → Fracture energy : 1.4868e+00 J
  → Total energy    : 3.2827e+00 J


## Step 121/796: t = 5.43e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 120 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.130e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.249e-04
  |ΔD|_∞ = 6.577e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.497e-05
  |ΔD|_∞ = 1.315e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.994e-06
  |ΔD|_∞ = 2.631e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.988e-07
  |ΔD|_∞ = 5.262e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.998e-07
  |ΔD|_∞ = 1.052e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.148e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.995e-08
  |ΔD|_∞ = 2.105e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7997e+00 J
  → Fracture energy : 1.4869e+00 J
  → Total energy    : 3.2867e+00 J


## Step 122/796: t = 5.48e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 121 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.01e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.129e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.251e-04
  |ΔD|_∞ = 6.605e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.01e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.502e-05
  |ΔD|_∞ = 1.321e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.01e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.495e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.004e-06
  |ΔD|_∞ = 2.642e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.01e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.495e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.001e-06
  |ΔD|_∞ = 5.284e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.01e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.989e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.001e-07
  |ΔD|_∞ = 1.057e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.01e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.989e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.003e-08
  |ΔD|_∞ = 2.113e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8036e+00 J
  → Fracture energy : 1.4871e+00 J
  → Total energy    : 3.2906e+00 J


## Step 123/796: t = 5.52e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 122 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.02e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.127e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.253e-04
  |ΔD|_∞ = 6.632e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.02e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.507e-05
  |ΔD|_∞ = 1.326e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.02e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.013e-06
  |ΔD|_∞ = 2.653e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.02e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.003e-06
  |ΔD|_∞ = 5.306e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.02e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.005e-07
  |ΔD|_∞ = 1.061e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.02e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.995e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.011e-08
  |ΔD|_∞ = 2.122e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8074e+00 J
  → Fracture energy : 1.4872e+00 J
  → Total energy    : 3.2946e+00 J


## Step 124/796: t = 5.57e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 123 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.03e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.126e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.256e-04
  |ΔD|_∞ = 6.659e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.03e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.043e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.511e-05
  |ΔD|_∞ = 1.332e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.03e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.022e-06
  |ΔD|_∞ = 2.664e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.03e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.004e-06
  |ΔD|_∞ = 5.328e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.03e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.009e-07
  |ΔD|_∞ = 1.066e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.03e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.018e-08
  |ΔD|_∞ = 2.131e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8113e+00 J
  → Fracture energy : 1.4873e+00 J
  → Total energy    : 3.2986e+00 J


## Step 125/796: t = 5.62e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 124 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.04e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.125e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.258e-04
  |ΔD|_∞ = 6.687e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.04e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.516e-05
  |ΔD|_∞ = 1.337e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.04e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.099e-20

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.032e-06
  |ΔD|_∞ = 2.675e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.04e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.099e-20

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.006e-06
  |ΔD|_∞ = 5.350e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.04e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.013e-07
  |ΔD|_∞ = 1.070e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.04e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.026e-08
  |ΔD|_∞ = 2.140e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8152e+00 J
  → Fracture energy : 1.4874e+00 J
  → Total energy    : 3.3026e+00 J


## Step 126/796: t = 5.66e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 125 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.05e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.124e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.260e-04
  |ΔD|_∞ = 6.716e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.05e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.521e-05
  |ΔD|_∞ = 1.343e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.05e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.042e-06
  |ΔD|_∞ = 2.686e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.05e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.246e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.008e-06
  |ΔD|_∞ = 5.373e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.05e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.716e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.017e-07
  |ΔD|_∞ = 1.075e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.05e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.203e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.034e-08
  |ΔD|_∞ = 2.149e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8191e+00 J
  → Fracture energy : 1.4875e+00 J
  → Total energy    : 3.3066e+00 J


## Step 127/796: t = 5.71e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 126 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.06e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.122e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.263e-04
  |ΔD|_∞ = 6.745e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.06e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.492e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.527e-05
  |ΔD|_∞ = 1.349e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.06e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.054e-06
  |ΔD|_∞ = 2.698e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.06e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.011e-06
  |ΔD|_∞ = 5.396e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.06e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.022e-07
  |ΔD|_∞ = 1.079e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.06e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.043e-08
  |ΔD|_∞ = 2.158e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8230e+00 J
  → Fracture energy : 1.4877e+00 J
  → Total energy    : 3.3106e+00 J


## Step 128/796: t = 5.75e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 127 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.07e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.121e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.266e-04
  |ΔD|_∞ = 6.774e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.07e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.490e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.531e-05
  |ΔD|_∞ = 1.355e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.07e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.047e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.062e-06
  |ΔD|_∞ = 2.709e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.07e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.047e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.012e-06
  |ΔD|_∞ = 5.419e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.07e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.025e-07
  |ΔD|_∞ = 1.084e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.07e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.029e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.050e-08
  |ΔD|_∞ = 2.168e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8268e+00 J
  → Fracture energy : 1.4878e+00 J
  → Total energy    : 3.3146e+00 J


## Step 129/796: t = 5.80e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 128 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.08e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.120e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.268e-04
  |ΔD|_∞ = 6.803e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.08e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.536e-05
  |ΔD|_∞ = 1.361e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.08e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.073e-06
  |ΔD|_∞ = 2.721e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.08e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.015e-06
  |ΔD|_∞ = 5.442e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.08e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.291e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.029e-07
  |ΔD|_∞ = 1.088e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.08e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.127e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.058e-08
  |ΔD|_∞ = 2.177e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8307e+00 J
  → Fracture energy : 1.4879e+00 J
  → Total energy    : 3.3186e+00 J


## Step 130/796: t = 5.84e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 129 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.09e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.119e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.271e-04
  |ΔD|_∞ = 6.833e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.09e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.272e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.543e-05
  |ΔD|_∞ = 1.367e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.09e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.085e-06
  |ΔD|_∞ = 2.733e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.09e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.212e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.017e-06
  |ΔD|_∞ = 5.466e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.09e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.212e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.034e-07
  |ΔD|_∞ = 1.093e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.09e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.068e-08
  |ΔD|_∞ = 2.187e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8346e+00 J
  → Fracture energy : 1.4880e+00 J
  → Total energy    : 3.3227e+00 J


## Step 131/796: t = 5.89e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 130 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.1e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.117e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.273e-04
  |ΔD|_∞ = 6.863e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.1e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.941e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.547e-05
  |ΔD|_∞ = 1.373e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.1e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.094e-06
  |ΔD|_∞ = 2.745e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.1e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.019e-06
  |ΔD|_∞ = 5.490e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.1e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.038e-07
  |ΔD|_∞ = 1.098e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.1e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.941e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.075e-08
  |ΔD|_∞ = 2.196e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8385e+00 J
  → Fracture energy : 1.4881e+00 J
  → Total energy    : 3.3267e+00 J


## Step 132/796: t = 5.93e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 131 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.11e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.116e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.276e-04
  |ΔD|_∞ = 6.893e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.11e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.552e-05
  |ΔD|_∞ = 1.379e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.11e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.105e-06
  |ΔD|_∞ = 2.757e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.11e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.021e-06
  |ΔD|_∞ = 5.514e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.11e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.042e-07
  |ΔD|_∞ = 1.103e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.11e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.084e-08
  |ΔD|_∞ = 2.206e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8424e+00 J
  → Fracture energy : 1.4883e+00 J
  → Total energy    : 3.3307e+00 J


## Step 133/796: t = 5.98e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 132 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.12e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.115e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.281e-04
  |ΔD|_∞ = 6.924e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.12e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.561e-05
  |ΔD|_∞ = 1.385e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.12e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.123e-06
  |ΔD|_∞ = 2.769e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.12e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.025e-06
  |ΔD|_∞ = 5.539e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.12e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.049e-07
  |ΔD|_∞ = 1.108e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.12e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.098e-08
  |ΔD|_∞ = 2.216e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8464e+00 J
  → Fracture energy : 1.4884e+00 J
  → Total energy    : 3.3347e+00 J


## Step 134/796: t = 6.02e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 133 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.13e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.114e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.282e-04
  |ΔD|_∞ = 6.954e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.13e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.563e-05
  |ΔD|_∞ = 1.391e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.13e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.126e-06
  |ΔD|_∞ = 2.782e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.13e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.025e-06
  |ΔD|_∞ = 5.563e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.13e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.051e-07
  |ΔD|_∞ = 1.113e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.13e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.101e-08
  |ΔD|_∞ = 2.225e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8503e+00 J
  → Fracture energy : 1.4885e+00 J
  → Total energy    : 3.3388e+00 J


## Step 135/796: t = 6.07e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 134 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.14e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.112e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.284e-04
  |ΔD|_∞ = 6.985e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.14e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.569e-05
  |ΔD|_∞ = 1.397e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.14e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.138e-06
  |ΔD|_∞ = 2.794e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.14e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.028e-06
  |ΔD|_∞ = 5.588e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.14e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.519e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.055e-07
  |ΔD|_∞ = 1.118e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.14e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.519e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.110e-08
  |ΔD|_∞ = 2.235e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8542e+00 J
  → Fracture energy : 1.4886e+00 J
  → Total energy    : 3.3428e+00 J


## Step 136/796: t = 6.11e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 135 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.15e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.111e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.287e-04
  |ΔD|_∞ = 7.017e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.15e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.574e-05
  |ΔD|_∞ = 1.403e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.15e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.149e-06
  |ΔD|_∞ = 2.807e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.15e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.030e-06
  |ΔD|_∞ = 5.613e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.15e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.060e-07
  |ΔD|_∞ = 1.123e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.15e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.119e-08
  |ΔD|_∞ = 2.245e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8581e+00 J
  → Fracture energy : 1.4887e+00 J
  → Total energy    : 3.3468e+00 J


## Step 137/796: t = 6.16e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 136 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.16e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.110e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.290e-04
  |ΔD|_∞ = 7.049e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.16e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.630e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.580e-05
  |ΔD|_∞ = 1.410e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.16e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.115e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.160e-06
  |ΔD|_∞ = 2.819e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.16e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.613e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.032e-06
  |ΔD|_∞ = 5.639e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.16e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.282e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.064e-07
  |ΔD|_∞ = 1.128e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.16e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.282e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.128e-08
  |ΔD|_∞ = 2.256e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8620e+00 J
  → Fracture energy : 1.4889e+00 J
  → Total energy    : 3.3509e+00 J


## Step 138/796: t = 6.20e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 137 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.17e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.109e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.293e-04
  |ΔD|_∞ = 7.081e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.17e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.114e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.586e-05
  |ΔD|_∞ = 1.416e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.17e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.172e-06
  |ΔD|_∞ = 2.832e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.17e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.115e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.034e-06
  |ΔD|_∞ = 5.665e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.17e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.115e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.069e-07
  |ΔD|_∞ = 1.133e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.17e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.114e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.138e-08
  |ΔD|_∞ = 2.266e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8659e+00 J
  → Fracture energy : 1.4890e+00 J
  → Total energy    : 3.3549e+00 J


## Step 139/796: t = 6.25e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 138 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.18e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.108e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.296e-04
  |ΔD|_∞ = 7.113e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.18e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.234e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.592e-05
  |ΔD|_∞ = 1.423e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.18e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.537e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.184e-06
  |ΔD|_∞ = 2.845e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.18e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.896e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.037e-06
  |ΔD|_∞ = 5.691e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.18e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.896e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.074e-07
  |ΔD|_∞ = 1.138e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.18e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.147e-08
  |ΔD|_∞ = 2.276e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8699e+00 J
  → Fracture energy : 1.4891e+00 J
  → Total energy    : 3.3590e+00 J


## Step 140/796: t = 6.29e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 139 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.19e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.106e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.299e-04
  |ΔD|_∞ = 7.146e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.19e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.772e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.598e-05
  |ΔD|_∞ = 1.429e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.19e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.196e-06
  |ΔD|_∞ = 2.858e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.19e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.039e-06
  |ΔD|_∞ = 5.717e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.19e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.979e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.078e-07
  |ΔD|_∞ = 1.143e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.19e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.005e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.157e-08
  |ΔD|_∞ = 2.287e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8738e+00 J
  → Fracture energy : 1.4892e+00 J
  → Total energy    : 3.3630e+00 J


## Step 141/796: t = 6.34e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 140 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.105e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.302e-04
  |ΔD|_∞ = 7.179e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.477e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.604e-05
  |ΔD|_∞ = 1.436e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.208e-06
  |ΔD|_∞ = 2.872e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.042e-06
  |ΔD|_∞ = 5.743e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.083e-07
  |ΔD|_∞ = 1.149e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.181e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.166e-08
  |ΔD|_∞ = 2.297e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8778e+00 J
  → Fracture energy : 1.4894e+00 J
  → Total energy    : 3.3671e+00 J


## Step 142/796: t = 6.38e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 141 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.21e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.104e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.305e-04
  |ΔD|_∞ = 7.213e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.21e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.038e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.610e-05
  |ΔD|_∞ = 1.443e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.21e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.220e-06
  |ΔD|_∞ = 2.885e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.21e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.044e-06
  |ΔD|_∞ = 5.770e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.21e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.088e-07
  |ΔD|_∞ = 1.154e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.21e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.176e-08
  |ΔD|_∞ = 2.308e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8817e+00 J
  → Fracture energy : 1.4895e+00 J
  → Total energy    : 3.3712e+00 J


## Step 143/796: t = 6.43e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 142 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.22e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.103e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.308e-04
  |ΔD|_∞ = 7.247e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.22e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.522e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.616e-05
  |ΔD|_∞ = 1.449e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.22e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.054e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.232e-06
  |ΔD|_∞ = 2.899e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.22e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.632e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.046e-06
  |ΔD|_∞ = 5.797e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.22e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.654e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.093e-07
  |ΔD|_∞ = 1.159e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.22e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.150e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.186e-08
  |ΔD|_∞ = 2.319e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8856e+00 J
  → Fracture energy : 1.4896e+00 J
  → Total energy    : 3.3752e+00 J


## Step 144/796: t = 6.48e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 143 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.23e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.101e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.311e-04
  |ΔD|_∞ = 7.281e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.23e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.794e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.622e-05
  |ΔD|_∞ = 1.456e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.23e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.894e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.245e-06
  |ΔD|_∞ = 2.912e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.23e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.184e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.049e-06
  |ΔD|_∞ = 5.825e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.23e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.205e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.098e-07
  |ΔD|_∞ = 1.165e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.23e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.133e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.196e-08
  |ΔD|_∞ = 2.330e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8896e+00 J
  → Fracture energy : 1.4897e+00 J
  → Total energy    : 3.3793e+00 J


## Step 145/796: t = 6.52e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 144 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.24e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.100e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.314e-04
  |ΔD|_∞ = 7.315e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.24e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.248e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.629e-05
  |ΔD|_∞ = 1.463e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.24e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.135e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.257e-06
  |ΔD|_∞ = 2.926e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.24e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.287e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.051e-06
  |ΔD|_∞ = 5.852e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.24e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.365e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.103e-07
  |ΔD|_∞ = 1.170e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.24e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.922e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.206e-08
  |ΔD|_∞ = 2.341e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8935e+00 J
  → Fracture energy : 1.4898e+00 J
  → Total energy    : 3.3834e+00 J


## Step 146/796: t = 6.57e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 145 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.25e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.099e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.317e-04
  |ΔD|_∞ = 7.350e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.25e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.336e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.635e-05
  |ΔD|_∞ = 1.470e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.25e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.654e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.270e-06
  |ΔD|_∞ = 2.940e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.25e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.822e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.054e-06
  |ΔD|_∞ = 5.880e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.25e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.647e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.108e-07
  |ΔD|_∞ = 1.176e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.25e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.216e-08
  |ΔD|_∞ = 2.352e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8975e+00 J
  → Fracture energy : 1.4900e+00 J
  → Total energy    : 3.3875e+00 J


## Step 147/796: t = 6.61e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 146 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.26e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.098e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.321e-04
  |ΔD|_∞ = 7.386e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.26e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.721e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.641e-05
  |ΔD|_∞ = 1.477e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.26e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.702e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.283e-06
  |ΔD|_∞ = 2.954e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.26e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.057e-06
  |ΔD|_∞ = 5.908e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.26e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.113e-07
  |ΔD|_∞ = 1.182e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.26e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.226e-08
  |ΔD|_∞ = 2.363e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9015e+00 J
  → Fracture energy : 1.4901e+00 J
  → Total energy    : 3.3916e+00 J


## Step 148/796: t = 6.66e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 147 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.27e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.097e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.324e-04
  |ΔD|_∞ = 7.421e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.27e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.270e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.648e-05
  |ΔD|_∞ = 1.484e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.27e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.296e-06
  |ΔD|_∞ = 2.969e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.27e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.270e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.059e-06
  |ΔD|_∞ = 5.937e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.27e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.270e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.118e-07
  |ΔD|_∞ = 1.187e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.27e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.237e-08
  |ΔD|_∞ = 2.375e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9054e+00 J
  → Fracture energy : 1.4902e+00 J
  → Total energy    : 3.3957e+00 J


## Step 149/796: t = 6.70e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 148 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.28e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.095e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.327e-04
  |ΔD|_∞ = 7.457e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.28e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.056e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.655e-05
  |ΔD|_∞ = 1.491e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.28e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.309e-06
  |ΔD|_∞ = 2.983e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.28e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.049e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.062e-06
  |ΔD|_∞ = 5.966e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.28e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.990e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.124e-07
  |ΔD|_∞ = 1.193e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.28e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.210e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.247e-08
  |ΔD|_∞ = 2.386e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9094e+00 J
  → Fracture energy : 1.4903e+00 J
  → Total energy    : 3.3997e+00 J


## Step 150/796: t = 6.75e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 149 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.29e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.094e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.331e-04
  |ΔD|_∞ = 7.494e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.29e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.175e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.661e-05
  |ΔD|_∞ = 1.499e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.29e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.196e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.322e-06
  |ΔD|_∞ = 2.998e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.29e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.196e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.064e-06
  |ΔD|_∞ = 5.995e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.29e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.218e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.129e-07
  |ΔD|_∞ = 1.199e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.29e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.218e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.258e-08
  |ΔD|_∞ = 2.398e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9134e+00 J
  → Fracture energy : 1.4905e+00 J
  → Total energy    : 3.4038e+00 J


## Step 151/796: t = 6.79e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 150 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.3e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.093e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.334e-04
  |ΔD|_∞ = 7.531e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.3e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.668e-05
  |ΔD|_∞ = 1.506e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.3e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.998e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.336e-06
  |ΔD|_∞ = 3.012e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.3e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.998e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.067e-06
  |ΔD|_∞ = 6.025e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.3e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.134e-07
  |ΔD|_∞ = 1.205e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.3e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.269e-08
  |ΔD|_∞ = 2.410e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9173e+00 J
  → Fracture energy : 1.4906e+00 J
  → Total energy    : 3.4079e+00 J


## Step 152/796: t = 6.84e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 151 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.31e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.092e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.337e-04
  |ΔD|_∞ = 7.568e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.31e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.364e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.675e-05
  |ΔD|_∞ = 1.514e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.31e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.693e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.350e-06
  |ΔD|_∞ = 3.027e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.31e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.226e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.070e-06
  |ΔD|_∞ = 6.054e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.31e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.291e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.140e-07
  |ΔD|_∞ = 1.211e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.31e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.291e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.280e-08
  |ΔD|_∞ = 2.422e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9213e+00 J
  → Fracture energy : 1.4907e+00 J
  → Total energy    : 3.4121e+00 J


## Step 153/796: t = 6.88e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 152 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.32e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.091e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.341e-04
  |ΔD|_∞ = 7.605e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.32e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.145e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.682e-05
  |ΔD|_∞ = 1.521e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.32e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.031e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.363e-06
  |ΔD|_∞ = 3.042e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.32e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.073e-06
  |ΔD|_∞ = 6.084e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.32e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.031e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.145e-07
  |ΔD|_∞ = 1.217e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.32e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.291e-08
  |ΔD|_∞ = 2.434e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9253e+00 J
  → Fracture energy : 1.4909e+00 J
  → Total energy    : 3.4162e+00 J


## Step 154/796: t = 6.93e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 153 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.33e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.090e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.344e-04
  |ΔD|_∞ = 7.644e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.33e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.387e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.689e-05
  |ΔD|_∞ = 1.529e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.33e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.387e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.377e-06
  |ΔD|_∞ = 3.057e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.33e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.075e-06
  |ΔD|_∞ = 6.115e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.33e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.439e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.151e-07
  |ΔD|_∞ = 1.223e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.33e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.439e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.302e-08
  |ΔD|_∞ = 2.446e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9293e+00 J
  → Fracture energy : 1.4910e+00 J
  → Total energy    : 3.4203e+00 J


## Step 155/796: t = 6.97e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 154 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.34e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.088e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.348e-04
  |ΔD|_∞ = 7.682e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.34e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.105e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.696e-05
  |ΔD|_∞ = 1.536e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.34e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.105e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.392e-06
  |ΔD|_∞ = 3.073e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.34e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.105e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.078e-06
  |ΔD|_∞ = 6.146e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.34e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.105e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.157e-07
  |ΔD|_∞ = 1.229e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.34e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.105e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.313e-08
  |ΔD|_∞ = 2.458e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9333e+00 J
  → Fracture energy : 1.4911e+00 J
  → Total energy    : 3.4244e+00 J


## Step 156/796: t = 7.02e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 155 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.35e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.087e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.351e-04
  |ΔD|_∞ = 7.721e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.35e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.072e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.703e-05
  |ΔD|_∞ = 1.544e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.35e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.406e-06
  |ΔD|_∞ = 3.088e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.35e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.081e-06
  |ΔD|_∞ = 6.177e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.35e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.162e-07
  |ΔD|_∞ = 1.235e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.35e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.938e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.325e-08
  |ΔD|_∞ = 2.471e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9373e+00 J
  → Fracture energy : 1.4912e+00 J
  → Total energy    : 3.4285e+00 J


## Step 157/796: t = 7.06e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 156 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.36e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.086e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.355e-04
  |ΔD|_∞ = 7.760e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.36e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.149e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.710e-05
  |ΔD|_∞ = 1.552e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.36e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.149e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.420e-06
  |ΔD|_∞ = 3.104e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.36e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.149e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.084e-06
  |ΔD|_∞ = 6.208e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.36e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.149e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.168e-07
  |ΔD|_∞ = 1.242e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.36e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.149e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.336e-08
  |ΔD|_∞ = 2.483e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9413e+00 J
  → Fracture energy : 1.4914e+00 J
  → Total energy    : 3.4326e+00 J


## Step 158/796: t = 7.11e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 157 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.37e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.085e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.359e-04
  |ΔD|_∞ = 7.800e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.37e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.406e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.718e-05
  |ΔD|_∞ = 1.560e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.37e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.421e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.435e-06
  |ΔD|_∞ = 3.120e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.37e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.421e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.087e-06
  |ΔD|_∞ = 6.240e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.37e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.174e-07
  |ΔD|_∞ = 1.248e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.37e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.348e-08
  |ΔD|_∞ = 2.496e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9453e+00 J
  → Fracture energy : 1.4915e+00 J
  → Total energy    : 3.4368e+00 J


## Step 159/796: t = 7.15e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 158 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.38e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.084e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.362e-04
  |ΔD|_∞ = 7.840e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.38e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.725e-05
  |ΔD|_∞ = 1.568e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.38e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.450e-06
  |ΔD|_∞ = 3.136e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.38e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.090e-06
  |ΔD|_∞ = 6.272e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.38e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.180e-07
  |ΔD|_∞ = 1.254e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.38e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.360e-08
  |ΔD|_∞ = 2.509e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9493e+00 J
  → Fracture energy : 1.4916e+00 J
  → Total energy    : 3.4409e+00 J


## Step 160/796: t = 7.20e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 159 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.39e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.083e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.366e-04
  |ΔD|_∞ = 7.881e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.39e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.732e-05
  |ΔD|_∞ = 1.576e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.39e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.465e-06
  |ΔD|_∞ = 3.152e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.39e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.093e-06
  |ΔD|_∞ = 6.305e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.39e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.186e-07
  |ΔD|_∞ = 1.261e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.39e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.372e-08
  |ΔD|_∞ = 2.522e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9533e+00 J
  → Fracture energy : 1.4917e+00 J
  → Total energy    : 3.4451e+00 J


## Step 161/796: t = 7.25e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 160 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.081e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.370e-04
  |ΔD|_∞ = 7.922e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.740e-05
  |ΔD|_∞ = 1.584e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.480e-06
  |ΔD|_∞ = 3.169e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.096e-06
  |ΔD|_∞ = 6.338e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.192e-07
  |ΔD|_∞ = 1.268e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.384e-08
  |ΔD|_∞ = 2.535e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9573e+00 J
  → Fracture energy : 1.4919e+00 J
  → Total energy    : 3.4492e+00 J


## Step 162/796: t = 7.29e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 161 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.41e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.080e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.374e-04
  |ΔD|_∞ = 7.963e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.41e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.747e-05
  |ΔD|_∞ = 1.593e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.41e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.495e-06
  |ΔD|_∞ = 3.185e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.41e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.099e-06
  |ΔD|_∞ = 6.371e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.41e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.198e-07
  |ΔD|_∞ = 1.274e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.41e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.396e-08
  |ΔD|_∞ = 2.548e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9613e+00 J
  → Fracture energy : 1.4920e+00 J
  → Total energy    : 3.4533e+00 J


## Step 163/796: t = 7.34e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 162 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.42e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.079e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.378e-04
  |ΔD|_∞ = 8.006e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.42e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.755e-05
  |ΔD|_∞ = 1.601e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.42e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.510e-06
  |ΔD|_∞ = 3.202e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.42e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.102e-06
  |ΔD|_∞ = 6.404e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.42e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.204e-07
  |ΔD|_∞ = 1.281e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.42e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.408e-08
  |ΔD|_∞ = 2.562e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9654e+00 J
  → Fracture energy : 1.4921e+00 J
  → Total energy    : 3.4575e+00 J


## Step 164/796: t = 7.38e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 163 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.43e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.078e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.382e-04
  |ΔD|_∞ = 8.048e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.43e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.763e-05
  |ΔD|_∞ = 1.610e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.43e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.526e-06
  |ΔD|_∞ = 3.219e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.43e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.105e-06
  |ΔD|_∞ = 6.438e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.43e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.210e-07
  |ΔD|_∞ = 1.288e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.43e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.421e-08
  |ΔD|_∞ = 2.575e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9694e+00 J
  → Fracture energy : 1.4923e+00 J
  → Total energy    : 3.4617e+00 J


## Step 165/796: t = 7.43e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 164 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.44e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.077e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.386e-04
  |ΔD|_∞ = 8.091e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.44e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.625e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.771e-05
  |ΔD|_∞ = 1.618e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.44e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.542e-06
  |ΔD|_∞ = 3.236e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.44e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.010e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.108e-06
  |ΔD|_∞ = 6.473e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.44e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.010e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.217e-07
  |ΔD|_∞ = 1.295e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.44e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.434e-08
  |ΔD|_∞ = 2.589e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9734e+00 J
  → Fracture energy : 1.4924e+00 J
  → Total energy    : 3.4658e+00 J


## Step 166/796: t = 7.47e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 165 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.45e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.076e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.390e-04
  |ΔD|_∞ = 8.135e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.45e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.779e-05
  |ΔD|_∞ = 1.627e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.45e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.558e-06
  |ΔD|_∞ = 3.254e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.45e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.634e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.112e-06
  |ΔD|_∞ = 6.508e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.45e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.634e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.223e-07
  |ΔD|_∞ = 1.302e-06

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.45e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.447e-08
  |ΔD|_∞ = 2.603e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9774e+00 J
  → Fracture energy : 1.4925e+00 J
  → Total energy    : 3.4700e+00 J


## Step 167/796: t = 7.52e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 166 | dt = 4.53e+00 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-07
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.46e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.075e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.394e-04
  |ΔD|_∞ = 8.179e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.46e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.787e-05
  |ΔD|_∞ = 1.636e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.46e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.574e-06
  |ΔD|_∞ = 3.271e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.46e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.115e-06
  |ΔD|_∞ = 6.543e-06

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.46e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
