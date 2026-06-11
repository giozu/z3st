Info    : Reading 'mesh.msh'...
Info    : 9 entities
Info    : 55 nodes
Info    : 68 elements
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
  → Time steps          : 51
  → Regime              : axisymmetric
  → Models active       :
      thermal    → ON
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
  → Temperature  : 1.0
  → Displacement : 1.0
  → Damage       : 0.4
  Adaptive relaxation disabled


[ThermalModel] initializer
[ThermalModel] options loaded from input.yaml:
  solver              : linear
  linear_solver       : direct_mumps
  rtol                : 1e-10
  stag_tol            : 1e-08
  convergence         : rel_norm
[MechanicalModel] initializer
[MechanicalModel] options loaded from input.yaml:
  solver              : nonlinear
  linear_solver       : direct_mumps
  rtol                : 1e-09
  stag_tol            : 1e-07
  convergence         : rel_norm
[spine.load_materials]
Material loaded: clad
  → k defined as constant: 17.0
  → Gc not defined for clad
  → constitutive model: lame
  → creep: Norton, A0 = 2.820e-24 Pa^-n/s, n = 3.00, Q = 1.200e+05 J/mol
[spine.initialize_fields]
[UPDATING q_third]

Initializing the temperature field...
  → Setting initial temperature for material: 'clad'
    Set 55 DOFs to 600.00 K
  Initial T: min=600.00 K, max=600.00 K, mean=600.00 K

Initializing the displacement field...
  Initial u: min=0.00e+00 m, max=0.00e+00 m, mean=0.00e+00 m



***


### spine - set_boundary_conditions


***



Loading boundary conditions from 'boundary_conditions.yaml'
  **[INFO]** Dirichlet thermal BC on 'clad' → 600.0 K (first step) at region 'outer'
  **[INFO]** Clamp_x mechanical BC on 'clad' → 0.0 (first step) at region 'inner'
  **[INFO]** Clamp_y mechanical BC on 'clad' → 0.0 (first step) at region 'bottom'
  **[INFO]** Dirichlet_y mechanical BC on 'clad' → 5e-05 (first step) at region 'top'
Computing symbolic result fields (strain, stress, ...)

**[INFO]** Hot-reload of allow-listed input.yaml parameters is active. Edit input.yaml during the run; changes apply at the next step boundary. Allowed keys: damage.{stag_tol,rtol,hybrid_constraint,gamma_star}, mechanical.{stag_tol,rtol}, thermal.{stag_tol,rtol}, solver_settings.{max_iters,relax_*}.


## Step 01/51: t = 0.00e+00 s | LHR = 0.00e+00 W/m

[UPDATING q_third]
  → dt=0: solving static step / initial condition
Computing symbolic result fields (strain, stress, ...)



***


### spine - solve


***



Current step = 0 | dt = 0.00e+00 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 7.062e-16

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.000e+00
  [creep] predictor rel change = 0.000e+00

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)


## Step 02/51: t = 2.00e+06 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 1 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.748e-03
  [creep] predictor rel change = 1.000e+00

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 7.790e-08
  [creep] predictor rel change = 2.835e-02

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 1.272e-06

Convergence check


#### Iteration 4/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
  [creep] clad: max equivalent creep strain = 1.3000e-04
Computing symbolic result fields (strain, stress, ...)


## Step 03/51: t = 4.00e+06 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 2 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.246e-03
  [creep] predictor rel change = 4.371e-01

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.764e-08
  [creep] predictor rel change = 2.425e-02

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 5.409e-07

Convergence check


#### Iteration 4/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
  [creep] clad: max equivalent creep strain = 2.2272e-04
Computing symbolic result fields (strain, stress, ...)


## Step 04/51: t = 6.00e+06 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 3 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 9.396e-04
  [creep] predictor rel change = 3.552e-01

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.135e-08
  [creep] predictor rel change = 2.112e-02

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 2.565e-07

Convergence check


#### Iteration 4/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
  [creep] clad: max equivalent creep strain = 2.9260e-04
Computing symbolic result fields (strain, stress, ...)


## Step 05/51: t = 8.00e+06 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 4 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 7.375e-04
  [creep] predictor rel change = 2.982e-01

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.215e-09
  [creep] predictor rel change = 1.866e-02

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 1.326e-07

Convergence check


#### Iteration 4/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
  [creep] clad: max equivalent creep strain = 3.4746e-04
Computing symbolic result fields (strain, stress, ...)


## Step 06/51: t = 1.00e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 5 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.969e-04
  [creep] predictor rel change = 2.564e-01

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.616e-09
  [creep] predictor rel change = 1.669e-02

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 7.346e-08

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 3.9186e-04
Computing symbolic result fields (strain, stress, ...)


## Step 07/51: t = 1.20e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 6 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 4.949e-04
  [creep] predictor rel change = 2.245e-01

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.408e-09
  [creep] predictor rel change = 1.508e-02

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 4.308e-08

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 4.2868e-04
Computing symbolic result fields (strain, stress, ...)


## Step 08/51: t = 1.40e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 7 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 4.184e-04
  [creep] predictor rel change = 1.995e-01

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 8.038e-10
  [creep] predictor rel change = 1.374e-02

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 2.650e-08

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 4.5980e-04
Computing symbolic result fields (strain, stress, ...)


## Step 09/51: t = 1.60e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 8 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.593e-04
  [creep] predictor rel change = 1.793e-01

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 4.815e-10
  [creep] predictor rel change = 1.261e-02

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 1.696e-08

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 4.8653e-04
Computing symbolic result fields (strain, stress, ...)


## Step 10/51: t = 1.80e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 9 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.126e-04
  [creep] predictor rel change = 1.627e-01

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.005e-10
  [creep] predictor rel change = 1.165e-02

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 1.123e-08

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 5.0979e-04
Computing symbolic result fields (strain, stress, ...)


## Step 11/51: t = 2.00e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 10 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.751e-04
  [creep] predictor rel change = 1.489e-01

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.943e-10
  [creep] predictor rel change = 1.081e-02

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 7.662e-09

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 5.3025e-04
Computing symbolic result fields (strain, stress, ...)


## Step 12/51: t = 2.20e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 11 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.444e-04
  [creep] predictor rel change = 1.372e-01

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.295e-10
  [creep] predictor rel change = 1.009e-02

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 5.362e-09

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 5.4843e-04
Computing symbolic result fields (strain, stress, ...)


## Step 13/51: t = 2.40e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 12 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.189e-04
  [creep] predictor rel change = 1.271e-01

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 8.858e-11
  [creep] predictor rel change = 9.457e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 3.838e-09

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 5.6471e-04
Computing symbolic result fields (strain, stress, ...)


## Step 14/51: t = 2.60e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 13 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.975e-04
  [creep] predictor rel change = 1.184e-01

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 6.205e-11
  [creep] predictor rel change = 8.895e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 2.803e-09

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 5.7940e-04
Computing symbolic result fields (strain, stress, ...)


## Step 15/51: t = 2.80e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 14 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.793e-04
  [creep] predictor rel change = 1.108e-01

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 4.437e-11
  [creep] predictor rel change = 8.395e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 2.083e-09

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 5.9274e-04
Computing symbolic result fields (strain, stress, ...)


## Step 16/51: t = 3.00e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 15 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.637e-04
  [creep] predictor rel change = 1.041e-01

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.232e-11
  [creep] predictor rel change = 7.947e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 1.573e-09

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 6.0491e-04
Computing symbolic result fields (strain, stress, ...)


## Step 17/51: t = 3.20e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 16 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.502e-04
  [creep] predictor rel change = 9.814e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.393e-11
  [creep] predictor rel change = 7.544e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 1.205e-09

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 6.1608e-04
Computing symbolic result fields (strain, stress, ...)


## Step 18/51: t = 3.40e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 17 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.384e-04
  [creep] predictor rel change = 9.282e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.799e-11
  [creep] predictor rel change = 7.178e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 9.353e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 6.2638e-04
Computing symbolic result fields (strain, stress, ...)


## Step 19/51: t = 3.60e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 18 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.281e-04
  [creep] predictor rel change = 8.803e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.371e-11
  [creep] predictor rel change = 6.846e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 7.344e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 6.3591e-04
Computing symbolic result fields (strain, stress, ...)


## Step 20/51: t = 3.80e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 19 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.190e-04
  [creep] predictor rel change = 8.371e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.058e-11
  [creep] predictor rel change = 6.543e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 5.829e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 6.4476e-04
Computing symbolic result fields (strain, stress, ...)


## Step 21/51: t = 4.00e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 20 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.109e-04
  [creep] predictor rel change = 7.978e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 8.254e-12
  [creep] predictor rel change = 6.265e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 4.673e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 6.5301e-04
Computing symbolic result fields (strain, stress, ...)


## Step 22/51: t = 4.20e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 21 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.036e-04
  [creep] predictor rel change = 7.621e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 6.507e-12
  [creep] predictor rel change = 6.009e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 3.780e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 6.6072e-04
Computing symbolic result fields (strain, stress, ...)


## Step 23/51: t = 4.40e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 22 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 9.716e-05
  [creep] predictor rel change = 7.293e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.180e-12
  [creep] predictor rel change = 5.773e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 3.083e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 6.6795e-04
Computing symbolic result fields (strain, stress, ...)


## Step 24/51: t = 4.60e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 23 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 9.131e-05
  [creep] predictor rel change = 6.992e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 4.159e-12
  [creep] predictor rel change = 5.554e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 2.535e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 6.7474e-04
Computing symbolic result fields (strain, stress, ...)


## Step 25/51: t = 4.80e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 24 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 8.603e-05
  [creep] predictor rel change = 6.715e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.368e-12
  [creep] predictor rel change = 5.352e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 2.098e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 6.8114e-04
Computing symbolic result fields (strain, stress, ...)


## Step 26/51: t = 5.00e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 25 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 8.123e-05
  [creep] predictor rel change = 6.458e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.747e-12
  [creep] predictor rel change = 5.163e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 1.749e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 6.8718e-04
Computing symbolic result fields (strain, stress, ...)


## Step 27/51: t = 5.20e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 26 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 7.685e-05
  [creep] predictor rel change = 6.220e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.257e-12
  [creep] predictor rel change = 4.987e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 1.467e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 6.9290e-04
Computing symbolic result fields (strain, stress, ...)


## Step 28/51: t = 5.40e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 27 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 7.286e-05
  [creep] predictor rel change = 5.999e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.867e-12
  [creep] predictor rel change = 4.822e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 1.238e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 6.9832e-04
Computing symbolic result fields (strain, stress, ...)


## Step 29/51: t = 5.60e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 28 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 6.919e-05
  [creep] predictor rel change = 5.793e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.553e-12
  [creep] predictor rel change = 4.668e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 1.050e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 7.0347e-04
Computing symbolic result fields (strain, stress, ...)


## Step 30/51: t = 5.80e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 29 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 6.582e-05
  [creep] predictor rel change = 5.600e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.300e-12
  [creep] predictor rel change = 4.524e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 8.949e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 7.0837e-04
Computing symbolic result fields (strain, stress, ...)


## Step 31/51: t = 6.00e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 30 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 6.271e-05
  [creep] predictor rel change = 5.420e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.094e-12
  [creep] predictor rel change = 4.387e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 7.667e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 7.1303e-04
Computing symbolic result fields (strain, stress, ...)


## Step 32/51: t = 6.20e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 31 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.984e-05
  [creep] predictor rel change = 5.250e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 9.256e-13
  [creep] predictor rel change = 4.259e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 6.598e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 7.1748e-04
Computing symbolic result fields (strain, stress, ...)


## Step 33/51: t = 6.40e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 32 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.717e-05
  [creep] predictor rel change = 5.091e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 7.867e-13
  [creep] predictor rel change = 4.138e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 5.703e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 7.2174e-04
Computing symbolic result fields (strain, stress, ...)


## Step 34/51: t = 6.60e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 33 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.470e-05
  [creep] predictor rel change = 4.941e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 6.718e-13
  [creep] predictor rel change = 4.024e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 4.949e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 7.2581e-04
Computing symbolic result fields (strain, stress, ...)


## Step 35/51: t = 6.80e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 34 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.240e-05
  [creep] predictor rel change = 4.800e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.761e-13
  [creep] predictor rel change = 3.915e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 4.311e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 7.2971e-04
Computing symbolic result fields (strain, stress, ...)


## Step 36/51: t = 7.00e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 35 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.026e-05
  [creep] predictor rel change = 4.667e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 4.961e-13
  [creep] predictor rel change = 3.813e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 3.769e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 7.3345e-04
Computing symbolic result fields (strain, stress, ...)


## Step 37/51: t = 7.20e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 36 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 4.825e-05
  [creep] predictor rel change = 4.540e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 4.289e-13
  [creep] predictor rel change = 3.715e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 3.306e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 7.3704e-04
Computing symbolic result fields (strain, stress, ...)


## Step 38/51: t = 7.40e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 37 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 4.638e-05
  [creep] predictor rel change = 4.420e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.721e-13
  [creep] predictor rel change = 3.622e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 2.911e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 7.4049e-04
Computing symbolic result fields (strain, stress, ...)


## Step 39/51: t = 7.60e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 38 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 4.462e-05
  [creep] predictor rel change = 4.307e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.240e-13
  [creep] predictor rel change = 3.534e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 2.570e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 7.4381e-04
Computing symbolic result fields (strain, stress, ...)


## Step 40/51: t = 7.80e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 39 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 4.297e-05
  [creep] predictor rel change = 4.199e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.830e-13
  [creep] predictor rel change = 3.450e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 2.276e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 7.4700e-04
Computing symbolic result fields (strain, stress, ...)


## Step 41/51: t = 8.00e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 40 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 4.142e-05
  [creep] predictor rel change = 4.096e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.480e-13
  [creep] predictor rel change = 3.370e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 2.021e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 7.5008e-04
Computing symbolic result fields (strain, stress, ...)


## Step 42/51: t = 8.20e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 41 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.996e-05
  [creep] predictor rel change = 3.998e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.180e-13
  [creep] predictor rel change = 3.293e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 1.800e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 7.5306e-04
Computing symbolic result fields (strain, stress, ...)


## Step 43/51: t = 8.40e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 42 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.858e-05
  [creep] predictor rel change = 3.905e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.922e-13
  [creep] predictor rel change = 3.220e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 1.607e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 7.5593e-04
Computing symbolic result fields (strain, stress, ...)


## Step 44/51: t = 8.60e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 43 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.728e-05
  [creep] predictor rel change = 3.815e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.699e-13
  [creep] predictor rel change = 3.150e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 1.438e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 7.5870e-04
Computing symbolic result fields (strain, stress, ...)


## Step 45/51: t = 8.80e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 44 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.605e-05
  [creep] predictor rel change = 3.730e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.506e-13
  [creep] predictor rel change = 3.083e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 1.290e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 7.6138e-04
Computing symbolic result fields (strain, stress, ...)


## Step 46/51: t = 9.00e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 45 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.489e-05
  [creep] predictor rel change = 3.649e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.339e-13
  [creep] predictor rel change = 3.019e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 1.160e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 7.6398e-04
Computing symbolic result fields (strain, stress, ...)


## Step 47/51: t = 9.20e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 46 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.378e-05
  [creep] predictor rel change = 3.571e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.193e-13
  [creep] predictor rel change = 2.957e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 1.045e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 7.6649e-04
Computing symbolic result fields (strain, stress, ...)


## Step 48/51: t = 9.40e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 47 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.274e-05
  [creep] predictor rel change = 3.496e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.065e-13
  [creep] predictor rel change = 2.898e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 9.440e-12

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 7.6893e-04
Computing symbolic result fields (strain, stress, ...)


## Step 49/51: t = 9.60e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 48 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.174e-05
  [creep] predictor rel change = 3.424e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 9.530e-14
  [creep] predictor rel change = 2.841e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 8.549e-12

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 7.7129e-04
Computing symbolic result fields (strain, stress, ...)


## Step 50/51: t = 9.80e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 49 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.080e-05
  [creep] predictor rel change = 3.355e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 8.547e-14
  [creep] predictor rel change = 2.786e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 7.744e-12

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 7.7358e-04
Computing symbolic result fields (strain, stress, ...)


## Step 51/51: t = 1.00e+08 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 50 | dt = 2.00e+06 s
Coupling = staggered
  → Max iterations              : 40
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.990e-05
  [creep] predictor rel change = 3.289e-02

Convergence check


#### Iteration 2/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 7.682e-14
  [creep] predictor rel change = 2.733e-03

Convergence check


#### Iteration 3/40


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 5e-05
  Building weak form, volume integrals (dx) for clad, tag = 10
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 7.041e-12

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
  [creep] clad: max equivalent creep strain = 7.7581e-04
Computing symbolic result fields (strain, stress, ...)

Simulation completed in 84.17 s
Total time steps solved: 51
