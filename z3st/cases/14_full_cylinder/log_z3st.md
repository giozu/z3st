Info    : Reading 'mesh.msh'...
Info    : 9 entities
Info    : 8181 nodes
Info    : 8360 elements
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
  → Time steps          : 1
  → Regime              : axisymmetric
  → Models active       :
      thermal    → ON
      mechanical → ON
      damage     → OFF
      cluster    → OFF
      plasticity → OFF
  → Gap conductance     : None (value = 0.0)



### FiniteElementSetup initializer

Mechanical element order: 1
Thermal function space (V_t): FunctionSpace(<Mesh #0>, Basix element (P, quadrilateral, 1, gll_warped, unset, False, float64, []))
Mechanical function space (V_m): FunctionSpace(<Mesh #0>, blocked element (Basix element (P, quadrilateral, 1, gll_warped, unset, False, float64, []), (2,)))
Scalar function space (Q): FunctionSpace(<Mesh #0>, Basix element (P, quadrilateral, 0, gll_warped, unset, True, float64, []))
[Solver] initializer
  Applied relaxation factor:
  → Temperature  : 0.8
  → Displacement : 0.6
  → Damage       : 0.4
  Adaptive relaxation disabled


[ThermalModel] initializer
[ThermalModel] options loaded from input.yaml:
  solver              : linear
  linear_solver       : iterative_hypre
  rtol                : 1e-06
  stag_tol            : 1e-06
  convergence         : rel_norm
[MechanicalModel] initializer
[MechanicalModel] options loaded from input.yaml:
  solver              : linear
  linear_solver       : iterative_hypre
  rtol                : 1e-06
  stag_tol            : 1e-06
  convergence         : rel_norm
[spine.load_materials]
Material loaded: oxide
  → k defined as symbolic function: materials.ceramic.k
  → Gc not defined for oxide
  → constitutive model: lame
  E               → 170000000000.0 (float)
  G               → 65891472868.21705 (float)
  T_ref           → 300.0 (float)
  _k_func         → <function k at 0x7b7decde3c40> (function)
  alpha           → 1.45e-05 (float)
  bulk_modulus    → 134920634920.6349 (float)
  constitutive_mode → lame (str)
  cp              → 330.0 (float)
  fissile         → True (bool)
  gamma_heating   → 0.0 (float)
  k               → materials.ceramic.k (str)
  lmbda           → 90992986341.82355 (float)
  mu_gamma        → 0.0 (float)
  name            → ceramic (str)
  nu              → 0.29 (float)
  rho             → 11000.0 (float)
[spine.initialize_fields]
[UPDATING q_third]
Fissile material
  q_third += 9.947e+04 W/m³ (fissile: True)
  Heat flux = 1.989e+03 W/m2

Initializing the temperature field...
  → Setting initial temperature for material: 'oxide'
    Set 8181 DOFs to 300.00 K
  Initial T: min=300.00 K, max=300.00 K, mean=300.00 K

Initializing the displacement field...
  Initial u: min=0.00e+00 m, max=0.00e+00 m, mean=0.00e+00 m

k expression for oxide → 2.5



***


### spine - set_boundary_conditions


***



Loading boundary conditions from 'boundary_conditions.yaml'
  **[INFO]** Dirichlet thermal BC on 'oxide' → 500.0 K at region 'outer'
  **[INFO]** Clamp_y mechanical BC on 'oxide' → 0.0 (first step) at region 'bottom'
Computing symbolic result fields (strain, stress, ...)


## Step 01/1: t = 0.00e+00 s | LHR = 5.00e+02 W/m

[UPDATING q_third]
Fissile material
  q_third += 9.947e+04 W/m³ (fissile: True)
  Heat flux = 1.989e+03 W/m2
  → dt=0: solving static step / initial condition
Computing symbolic result fields (strain, stress, ...)



***


### spine - solve


***



Current step = 0 | dt = 0.00e+00 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for oxide, tag = 10
  → q_third[oxide](W/m3) min = 9.95e+04, max = 9.95e+04, mean = 9.95e+04
  Linear solver
  T_new: min=500.00 K, max=515.92 K, mean=510.58 K
  ||ΔT||/||T|| = 3.603e-01

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for oxide, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.000e+00

Convergence check


#### Iteration 2/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for oxide, tag = 10
  → q_third[oxide](W/m3) min = 9.95e+04, max = 9.95e+04, mean = 9.95e+04
  Linear solver
  T_new: min=500.00 K, max=515.92 K, mean=510.58 K
  ||ΔT||/||T|| = 7.156e-02

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for oxide, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.982e-01

Convergence check


#### Iteration 3/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for oxide, tag = 10
  → q_third[oxide](W/m3) min = 9.95e+04, max = 9.95e+04, mean = 9.95e+04
  Linear solver
  T_new: min=500.00 K, max=515.92 K, mean=510.58 K
  ||ΔT||/||T|| = 1.431e-02

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for oxide, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.789e-01

Convergence check


#### Iteration 4/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for oxide, tag = 10
  → q_third[oxide](W/m3) min = 9.95e+04, max = 9.95e+04, mean = 9.95e+04
  Linear solver
  T_new: min=500.00 K, max=515.92 K, mean=510.58 K
  ||ΔT||/||T|| = 2.862e-03

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for oxide, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.195e-01

Convergence check


#### Iteration 5/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for oxide, tag = 10
  → q_third[oxide](W/m3) min = 9.95e+04, max = 9.95e+04, mean = 9.95e+04
  Linear solver
  T_new: min=500.00 K, max=515.92 K, mean=510.58 K
  ||ΔT||/||T|| = 5.724e-04

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for oxide, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.938e-02

Convergence check


#### Iteration 6/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for oxide, tag = 10
  → q_third[oxide](W/m3) min = 9.95e+04, max = 9.95e+04, mean = 9.95e+04
  Linear solver
  T_new: min=500.00 K, max=515.92 K, mean=510.58 K
  ||ΔT||/||T|| = 1.145e-04

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for oxide, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.007e-02

Convergence check


#### Iteration 7/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for oxide, tag = 10
  → q_third[oxide](W/m3) min = 9.95e+04, max = 9.95e+04, mean = 9.95e+04
  Linear solver
  T_new: min=500.00 K, max=515.92 K, mean=510.58 K
  ||ΔT||/||T|| = 2.290e-05

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for oxide, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.091e-03

Convergence check


#### Iteration 8/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for oxide, tag = 10
  → q_third[oxide](W/m3) min = 9.95e+04, max = 9.95e+04, mean = 9.95e+04
  Linear solver
  T_new: min=500.00 K, max=515.92 K, mean=510.58 K
  ||ΔT||/||T|| = 4.580e-06

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for oxide, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.249e-03

Convergence check


#### Iteration 9/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for oxide, tag = 10
  → q_third[oxide](W/m3) min = 9.95e+04, max = 9.95e+04, mean = 9.95e+04
  Linear solver
  T_new: min=500.00 K, max=515.92 K, mean=510.58 K
  ||ΔT||/||T|| = 9.159e-07

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for oxide, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.302e-03

Convergence check


#### Iteration 10/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for oxide, tag = 10
  → q_third[oxide](W/m3) min = 9.95e+04, max = 9.95e+04, mean = 9.95e+04
  Linear solver
  T_new: min=500.00 K, max=515.92 K, mean=510.58 K
  ||ΔT||/||T|| = 1.832e-07

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for oxide, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.214e-04

Convergence check


#### Iteration 11/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for oxide, tag = 10
  → q_third[oxide](W/m3) min = 9.95e+04, max = 9.95e+04, mean = 9.95e+04
  Linear solver
  T_new: min=500.00 K, max=515.92 K, mean=510.58 K
  ||ΔT||/||T|| = 3.664e-08

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for oxide, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.087e-04

Convergence check


#### Iteration 12/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for oxide, tag = 10
  → q_third[oxide](W/m3) min = 9.95e+04, max = 9.95e+04, mean = 9.95e+04
  Linear solver
  T_new: min=500.00 K, max=515.92 K, mean=510.58 K
  ||ΔT||/||T|| = 7.327e-09

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for oxide, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.349e-05

Convergence check


#### Iteration 13/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for oxide, tag = 10
  → q_third[oxide](W/m3) min = 9.95e+04, max = 9.95e+04, mean = 9.95e+04
  Linear solver
  T_new: min=500.00 K, max=515.92 K, mean=510.58 K
  ||ΔT||/||T|| = 1.465e-09

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for oxide, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.340e-05

Convergence check


#### Iteration 14/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for oxide, tag = 10
  → q_third[oxide](W/m3) min = 9.95e+04, max = 9.95e+04, mean = 9.95e+04
  Linear solver
  T_new: min=500.00 K, max=515.92 K, mean=510.58 K
  ||ΔT||/||T|| = 2.931e-10

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for oxide, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.336e-05

Convergence check


#### Iteration 15/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for oxide, tag = 10
  → q_third[oxide](W/m3) min = 9.95e+04, max = 9.95e+04, mean = 9.95e+04
  Linear solver
  T_new: min=500.00 K, max=515.92 K, mean=510.58 K
  ||ΔT||/||T|| = 5.862e-11

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for oxide, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.344e-06

Convergence check


#### Iteration 16/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for oxide, tag = 10
  → q_third[oxide](W/m3) min = 9.95e+04, max = 9.95e+04, mean = 9.95e+04
  Linear solver
  T_new: min=500.00 K, max=515.92 K, mean=510.58 K
  ||ΔT||/||T|| = 1.172e-11

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for oxide, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.138e-06

Convergence check


#### Iteration 17/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for oxide, tag = 10
  → q_third[oxide](W/m3) min = 9.95e+04, max = 9.95e+04, mean = 9.95e+04
  Linear solver
  T_new: min=500.00 K, max=515.92 K, mean=510.58 K
  ||ΔT||/||T|| = 2.345e-12

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for oxide, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.551e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 17 iterations.
Computing symbolic result fields (strain, stress, ...)

Simulation completed in 3.41 s
Total time steps solved: 1
