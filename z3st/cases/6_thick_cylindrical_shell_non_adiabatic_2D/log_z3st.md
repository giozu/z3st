Info    : Reading 'mesh.msh'...
Info    : 9 entities
Info    : 3321 nodes
Info    : 3440 elements
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
  → relax_shrink : 0.7
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
  debug               : False
[spine.load_materials]
Material loaded: steel
  → k defined as constant: 48.1
  → Gc not defined for steel
  → constitutive model: lame
  E               → 177000000000.0 (float)
  G               → 68076923076.92307 (float)
  T_ref           → 300.0 (float)
  alpha           → 1.7e-05 (float)
  bulk_modulus    → 147499999999.99997 (float)
  constitutive_mode → lame (str)
  cp              → 200.0 (float)
  gamma_heating   → 2000000.0 (float)
  k               → 48.1 (float)
  lmbda           → 102115384615.38461 (float)
  mu_gamma        → 24.0 (float)
  name            → vessel_steel (str)
  nu              → 0.3 (float)
  rho             → 8000.0 (float)
[spine.initialize_fields]
[UPDATING q_third]

Initializing the temperature field...
  → Setting initial temperature for material: 'steel'
    Set 3321 DOFs to 300.00 K
  Initial T: min=300.00 K, max=300.00 K, mean=300.00 K

Initializing the displacement field...
  Initial u: min=0.00e+00 m, max=0.00e+00 m, mean=0.00e+00 m



***


### spine - set_boundary_conditions


***



Loading boundary conditions from 'boundary_conditions.yaml'
  **[INFO]** Dirichlet thermal BC on 'steel' → 490.0 K at region 'inner_radius'
  **[INFO]** Dirichlet thermal BC on 'steel' → 500.0 K at region 'outer_radius'
  **[INFO]** Clamp_y mechanical BC on 'steel' → 0.0 (first step) at region 'bottom'


## Step 01/1: t = 0.00e+00 s | LHR = 0.00e+00 W/m

[UPDATING q_third]
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

  Building weak form, volume integrals (dx) for steel, tag = 10
  → q_third[steel](W/m3) min = 1.24e+02, max = 2.00e+06, mean = 2.27e+05
  Linear solver
  T_new: min=490.00 K, max=539.38 K, mean=522.13 K
  ||ΔT||/||T|| = 4.015e-01
  [adaptive] relax_T=0.90

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.000e+00
  [adaptive] relax_u=0.40

Convergence check


#### Iteration 2/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 10
  → q_third[steel](W/m3) min = 1.24e+02, max = 2.00e+06, mean = 2.27e+05
  Linear solver
  T_new: min=490.00 K, max=539.38 K, mean=522.13 K
  ||ΔT||/||T|| = 3.922e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.976e-01
  [adaptive] relax_u=0.44

Convergence check


#### Iteration 3/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 10
  → q_third[steel](W/m3) min = 1.24e+02, max = 2.00e+06, mean = 2.27e+05
  Linear solver
  T_new: min=490.00 K, max=539.38 K, mean=522.13 K
  ||ΔT||/||T|| = 4.140e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.717e-01
  [adaptive] relax_u=0.48

Convergence check


#### Iteration 4/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 10
  → q_third[steel](W/m3) min = 1.24e+02, max = 2.00e+06, mean = 2.27e+05
  Linear solver
  T_new: min=490.00 K, max=539.38 K, mean=522.13 K
  ||ΔT||/||T|| = 2.070e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.912e-01
  [adaptive] relax_u=0.53

Convergence check


#### Iteration 5/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 10
  → q_third[steel](W/m3) min = 1.24e+02, max = 2.00e+06, mean = 2.27e+05
  Linear solver
  T_new: min=490.00 K, max=539.38 K, mean=522.13 K
  ||ΔT||/||T|| = 1.035e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.653e-01
  [adaptive] relax_u=0.59

Convergence check


#### Iteration 6/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 10
  → q_third[steel](W/m3) min = 1.24e+02, max = 2.00e+06, mean = 2.27e+05
  Linear solver
  T_new: min=490.00 K, max=539.38 K, mean=522.13 K
  ||ΔT||/||T|| = 5.175e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.504e-02
  [adaptive] relax_u=0.64

Convergence check


#### Iteration 7/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 10
  → q_third[steel](W/m3) min = 1.24e+02, max = 2.00e+06, mean = 2.27e+05
  Linear solver
  T_new: min=490.00 K, max=539.38 K, mean=522.13 K
  ||ΔT||/||T|| = 2.587e-08
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.876e-02
  [adaptive] relax_u=0.71

Convergence check


#### Iteration 8/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 10
  → q_third[steel](W/m3) min = 1.24e+02, max = 2.00e+06, mean = 2.27e+05
  Linear solver
  T_new: min=490.00 K, max=539.38 K, mean=522.13 K
  ||ΔT||/||T|| = 1.294e-09
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.517e-02
  [adaptive] relax_u=0.78

Convergence check


#### Iteration 9/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 10
  → q_third[steel](W/m3) min = 1.24e+02, max = 2.00e+06, mean = 2.27e+05
  Linear solver
  T_new: min=490.00 K, max=539.38 K, mean=522.13 K
  ||ΔT||/||T|| = 6.469e-11
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.862e-03
  [adaptive] relax_u=0.86

Convergence check


#### Iteration 10/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 10
  → q_third[steel](W/m3) min = 1.24e+02, max = 2.00e+06, mean = 2.27e+05
  Linear solver
  T_new: min=490.00 K, max=539.38 K, mean=522.13 K
  ||ΔT||/||T|| = 3.234e-12
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.179e-03
  [adaptive] relax_u=0.94

Convergence check


#### Iteration 11/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 10
  → q_third[steel](W/m3) min = 1.24e+02, max = 2.00e+06, mean = 2.27e+05
  Linear solver
  T_new: min=490.00 K, max=539.38 K, mean=522.13 K
  ||ΔT||/||T|| = 1.617e-13
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.849e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 12/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 10
  → q_third[steel](W/m3) min = 1.24e+02, max = 2.00e+06, mean = 2.27e+05
  Linear solver
  T_new: min=490.00 K, max=539.38 K, mean=522.13 K
  ||ΔT||/||T|| = 8.084e-15
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.058e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 13/100


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for steel, tag = 10
  → q_third[steel](W/m3) min = 1.24e+02, max = 2.00e+06, mean = 2.27e+05
  Linear solver
  T_new: min=490.00 K, max=539.38 K, mean=522.13 K
  ||ΔT||/||T|| = 4.284e-16
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.292e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 13 iterations.
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
VTU file exported to: output/fields.vtu

Simulation completed in 2.88 s
Total time steps solved: 1
