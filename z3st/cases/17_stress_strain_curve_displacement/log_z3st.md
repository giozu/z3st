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
  → Time steps          : 5
  → Regime              : 2d
  → Models active       :
      thermal    → OFF
      mechanical → ON
      damage     → OFF
      cluster    → OFF
      plasticity → OFF
      contact    → OFF
  → Gap conductance     : None (value = 0.0)



### FiniteElementSetup initializer

Mechanical element order: 1
Thermal function space (V_t): FunctionSpace(<Mesh #0>, Basix element (P, quadrilateral, 1, gll_warped, unset, False, float64, []))
Mechanical function space (V_m): FunctionSpace(<Mesh #0>, blocked element (Basix element (P, quadrilateral, 1, gll_warped, unset, False, float64, []), (2,)))
Scalar function space (Q): FunctionSpace(<Mesh #0>, Basix element (P, quadrilateral, 0, gll_warped, unset, True, float64, []))
[Solver] initializer
  Applied relaxation factor:
  → Temperature  : 0.9
  → Displacement : 0.7
  → Damage       : 0.4
  Adaptive relaxation enabled
  → relax_growth  : 1.1
  → relax_shrink : 0.9
  → relax_min  : 0.05
  → relax_max : 0.95


[MechanicalModel] initializer
[MechanicalModel] options loaded from input.yaml:
  solver              : linear
  linear_solver       : iterative_hypre
  rtol                : 1e-08
  stag_tol            : 1e-06
  convergence         : rel_norm
  debug               : False
[spine.load_materials]
Material loaded: steel
  → k defined as constant: 50.0
  → Gc not defined for steel
  → constitutive model: lame
  E               → 200000000000.0 (float)
  G               → 76923076923.07692 (float)
  T_ref           → 300.0 (float)
  alpha           → 1e-05 (float)
  bulk_modulus    → 166666666666.66666 (float)
  constitutive_mode → lame (str)
  cp              → 200.0 (float)
  gamma_heating   → 0.0 (float)
  hardening_modulus → 10000000000.0 (float)
  k               → 50.0 (float)
  lmbda           → 115384615384.61539 (float)
  mu_gamma        → 25.0 (float)
  name            → steel (str)
  nu              → 0.3 (float)
  rho             → 8000.0 (float)
  sigma_c         → 600000000.0 (float)
  yield_strength  → 200000000.0 (float)
[spine.initialize_fields]

Initializing the displacement field...
  Initial u: min=0.00e+00 m, max=0.00e+00 m, mean=0.00e+00 m



***


### spine - set_boundary_conditions


***



Loading boundary conditions from 'boundary_conditions.yaml'
  **[INFO]** Clamp_x mechanical BC on 'steel' → 0.0 (first step) at region 'xmin'
  **[INFO]** Step-dependent Dirichlet list (2D), length 5
  **[INFO]** Dirichlet mechanical BC on 'steel' → [0.0, 0.0] at region 'xmax'
Computing symbolic result fields (strain, stress, ...)

**[INFO]** Hot-reload of allow-listed input.yaml parameters is active. Edit input.yaml during the run; changes apply at the next step boundary. Allowed keys: damage.{stag_tol,rtol,hybrid_constraint,gamma_star}, mechanical.{stag_tol,rtol}, thermal.{stag_tol,rtol}, solver_settings.{max_iters,relax_*}.


## Step 01/5: t = 0.00e+00 s | LHR = 0.00e+00 W/m

  → dt=0: solving static step / initial condition
Computing symbolic result fields (strain, stress, ...)



***


### spine - solve


***



Current step = 0 | dt = 0.00e+00 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-08
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [0.0, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 5
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=0.70

Convergence check

**[SUCCESS]** Staggered solver converged in 1 iterations.
Computing symbolic result fields (strain, stress, ...)


## Step 02/5: t = 2.50e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 1 | dt = 2.50e+00 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-08
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [5e-11, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 5
  Linear solver
  ||Δu||/||u|| = 1.000e+00
  [adaptive] relax_u=0.70

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [5e-11, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 5
  Linear solver
  ||Δu||/||u|| = 2.890e-01
  [adaptive] relax_u=0.77

Convergence check


#### Iteration 3/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [5e-11, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 5
  Linear solver
  ||Δu||/||u|| = 9.536e-02
  [adaptive] relax_u=0.85

Convergence check


#### Iteration 4/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [5e-11, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 5
  Linear solver
  ||Δu||/||u|| = 2.413e-02
  [adaptive] relax_u=0.93

Convergence check


#### Iteration 5/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [5e-11, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 5
  Linear solver
  ||Δu||/||u|| = 4.060e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 6/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [5e-11, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 5
  Linear solver
  ||Δu||/||u|| = 2.828e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 7/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [5e-11, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 5
  Linear solver
  ||Δu||/||u|| = 1.414e-05
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 8/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [5e-11, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 5
  Linear solver
  ||Δu||/||u|| = 7.044e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)


## Step 03/5: t = 5.00e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 2 | dt = 2.50e+00 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-08
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [5e-09, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 5
  Linear solver
  ||Δu||/||u|| = 9.895e-01
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [5e-09, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 5
  Linear solver
  ||Δu||/||u|| = 4.846e-02
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [5e-09, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 5
  Linear solver
  ||Δu||/||u|| = 2.423e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [5e-09, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 5
  Linear solver
  ||Δu||/||u|| = 1.212e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 5/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [5e-09, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 5
  Linear solver
  ||Δu||/||u|| = 6.058e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 6/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [5e-09, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 5
  Linear solver
  ||Δu||/||u|| = 3.029e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)


## Step 04/5: t = 7.50e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 3 | dt = 2.50e+00 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-08
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [5e-07, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 5
  Linear solver
  ||Δu||/||u|| = 9.895e-01
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [5e-07, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 5
  Linear solver
  ||Δu||/||u|| = 4.846e-02
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [5e-07, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 5
  Linear solver
  ||Δu||/||u|| = 2.423e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [5e-07, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 5
  Linear solver
  ||Δu||/||u|| = 1.212e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 5/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [5e-07, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 5
  Linear solver
  ||Δu||/||u|| = 6.058e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 6/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [5e-07, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 5
  Linear solver
  ||Δu||/||u|| = 3.029e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)


## Step 05/5: t = 1.00e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 4 | dt = 2.50e+00 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-08
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [5e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 5
  Linear solver
  ||Δu||/||u|| = 9.895e-01
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [5e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 5
  Linear solver
  ||Δu||/||u|| = 4.846e-02
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 3/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [5e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 5
  Linear solver
  ||Δu||/||u|| = 2.423e-03
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 4/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [5e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 5
  Linear solver
  ||Δu||/||u|| = 1.212e-04
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 5/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [5e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 5
  Linear solver
  ||Δu||/||u|| = 6.057e-06
  [adaptive] relax_u=0.95

Convergence check


#### Iteration 6/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 2 → [5e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 5
  Linear solver
  ||Δu||/||u|| = 3.026e-07
  [adaptive] relax_u=0.95

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)

Simulation completed in 33.07 s
Total time steps solved: 5
