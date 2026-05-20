Info    : Reading 'mesh.msh'...
Info    : 9 entities
Info    : 6561 nodes
Info    : 6720 elements
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
  → Temperature  : 0.9
  → Displacement : 0.4
  → Damage       : 0.4
  Adaptive relaxation enabled
  → relax_growth  : 1.1
  → relax_shrink : 0.9
  → relax_min  : 0.05
  → relax_max : 0.95


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
Material loaded: volume
  → k defined as symbolic function: materials.ceramic.k
  → Gc not defined for volume
  → constitutive model: lame
  E               → 170000000000.0 (float)
  G               → 65891472868.21705 (float)
  T_ref           → 300.0 (float)
  _k_func         → <function k at 0x70d4bbd17ce0> (function)
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
  q_third += 1.005e+05 W/m³ (fissile: True)
  Heat flux = 1.989e+03 W/m2

Initializing the temperature field...
  → Setting initial temperature for material: 'volume'
    Set 6561 DOFs to 300.00 K
  Initial T: min=300.00 K, max=300.00 K, mean=300.00 K

Initializing the displacement field...
  Initial u: min=0.00e+00 m, max=0.00e+00 m, mean=0.00e+00 m

k expression for volume → 2.5



***


### spine - set_boundary_conditions


***



Loading boundary conditions from 'boundary_conditions.yaml'
  **[INFO]** Dirichlet thermal BC on 'volume' → 500.0 K at region 'outer'
  **[INFO]** Clamp_y mechanical BC on 'volume' → 0.0 (first step) at region 'bottom'
Computing symbolic result fields (strain, stress, ...)

**[INFO]** Hot-reload of allow-listed input.yaml parameters is active. Edit input.yaml during the run; changes apply at the next step boundary. Allowed keys: damage.{stag_tol,rtol,hybrid_constraint,gamma_star}, mechanical.{stag_tol,rtol}, thermal.{stag_tol,rtol}, solver_settings.{max_iters,relax_*}.


## Step 01/1: t = 0.00e+00 s | LHR = 5.00e+02 W/m

[UPDATING q_third]
Fissile material
  q_third += 1.005e+05 W/m³ (fissile: True)
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

  Building weak form, volume integrals (dx) for volume, tag = 10
  → q_third[volume](W/m3) min = 1.00e+05, max = 1.00e+05, mean = 1.00e+05
  Linear solver
  T_new: min=500.00 K, max=515.18 K, mean=509.86 K
  ||ΔT||/||T|| = 3.867e-01
  [adaptive] relax_T=0.90

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for volume, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.000e+00
  [adaptive] relax_u=0.40

Convergence check


#### Iteration 2/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for volume, tag = 10
  → q_third[volume](W/m3) min = 1.00e+05, max = 1.00e+05, mean = 1.00e+05
  Linear solver
  T_new: min=500.00 K, max=515.18 K, mean=509.86 K
  ||ΔT||/||T|| = 3.841e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for volume, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.991e-01
  [adaptive] relax_u=0.44

Convergence check


#### Iteration 3/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for volume, tag = 10
  → q_third[volume](W/m3) min = 1.00e+05, max = 1.00e+05, mean = 1.00e+05
  Linear solver
  T_new: min=500.00 K, max=515.18 K, mean=509.86 K
  ||ΔT||/||T|| = 4.054e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for volume, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.729e-01
  [adaptive] relax_u=0.48

Convergence check


#### Iteration 4/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for volume, tag = 10
  → q_third[volume](W/m3) min = 1.00e+05, max = 1.00e+05, mean = 1.00e+05
  Linear solver
  T_new: min=500.00 K, max=515.18 K, mean=509.86 K
  ||ΔT||/||T|| = 2.027e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for volume, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.919e-01
  [adaptive] relax_u=0.53

Convergence check


#### Iteration 5/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for volume, tag = 10
  → q_third[volume](W/m3) min = 1.00e+05, max = 1.00e+05, mean = 1.00e+05
  Linear solver
  T_new: min=500.00 K, max=515.18 K, mean=509.86 K
  ||ΔT||/||T|| = 1.013e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for volume, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.657e-01
  [adaptive] relax_u=0.59

Convergence check


#### Iteration 6/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for volume, tag = 10
  → q_third[volume](W/m3) min = 1.00e+05, max = 1.00e+05, mean = 1.00e+05
  Linear solver
  T_new: min=500.00 K, max=515.18 K, mean=509.86 K
  ||ΔT||/||T|| = 5.067e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for volume, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.525e-02
  [adaptive] relax_u=0.64

Convergence check


#### Iteration 7/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for volume, tag = 10
  → q_third[volume](W/m3) min = 1.00e+05, max = 1.00e+05, mean = 1.00e+05
  Linear solver
  T_new: min=500.00 K, max=515.18 K, mean=509.86 K
  ||ΔT||/||T|| = 2.534e-08
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for volume, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.886e-02
  [adaptive] relax_u=0.71

Convergence check


#### Iteration 8/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for volume, tag = 10
  → q_third[volume](W/m3) min = 1.00e+05, max = 1.00e+05, mean = 1.00e+05
  Linear solver
  T_new: min=500.00 K, max=515.18 K, mean=509.86 K
  ||ΔT||/||T|| = 1.267e-09
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for volume, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.521e-02
  [adaptive] relax_u=0.78

Convergence check


#### Iteration 9/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for volume, tag = 10
  → q_third[volume](W/m3) min = 1.00e+05, max = 1.00e+05, mean = 1.00e+05
  Linear solver
  T_new: min=500.00 K, max=515.18 K, mean=509.86 K
  ||ΔT||/||T|| = 6.334e-11
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for volume, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.874e-03
  [adaptive] relax_u=0.86

Convergence check


#### Iteration 10/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for volume, tag = 10
  → q_third[volume](W/m3) min = 1.00e+05, max = 1.00e+05, mean = 1.00e+05
  Linear solver
  T_new: min=500.00 K, max=515.18 K, mean=509.86 K
  ||ΔT||/||T|| = 3.167e-12
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for volume, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.182e-03
  [adaptive] relax_u=0.94

Convergence check


#### Iteration 11/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for volume, tag = 10
  → q_third[volume](W/m3) min = 1.00e+05, max = 1.00e+05, mean = 1.00e+05
  Linear solver
  T_new: min=500.00 K, max=515.18 K, mean=509.86 K
  ||ΔT||/||T|| = 1.583e-13
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for volume, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.854e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 12/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for volume, tag = 10
  → q_third[volume](W/m3) min = 1.00e+05, max = 1.00e+05, mean = 1.00e+05
  Linear solver
  T_new: min=500.00 K, max=515.18 K, mean=509.86 K
  ||ΔT||/||T|| = 7.915e-15
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for volume, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.061e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 13/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for volume, tag = 10
  → q_third[volume](W/m3) min = 1.00e+05, max = 1.00e+05, mean = 1.00e+05
  Linear solver
  T_new: min=500.00 K, max=515.18 K, mean=509.86 K
  ||ΔT||/||T|| = 4.249e-16
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for volume, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.306e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 13 iterations.
Computing symbolic result fields (strain, stress, ...)

Simulation completed in 3.02 s
Total time steps solved: 1
