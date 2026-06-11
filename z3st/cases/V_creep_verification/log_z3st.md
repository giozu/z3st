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
  → Time steps          : 11
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
  **[INFO]** Neumann mechanical BC on 'clad' → top: 100000000.0 Pa (list loaded)
Computing symbolic result fields (strain, stress, ...)

**[INFO]** Hot-reload of allow-listed input.yaml parameters is active. Edit input.yaml during the run; changes apply at the next step boundary. Allowed keys: damage.{stag_tol,rtol,hybrid_constraint,gamma_star}, mechanical.{stag_tol,rtol}, thermal.{stag_tol,rtol}, solver_settings.{max_iters,relax_*}.


## Step 01/11: t = 0.00e+00 s | LHR = 0.00e+00 W/m

[UPDATING q_third]
  → dt=0: solving static step / initial condition
Computing symbolic result fields (strain, stress, ...)



***


### spine - solve


***



Current step = 0 | dt = 0.00e+00 s
Coupling = staggered
  → Max iterations              : 30
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 7.062e-16

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.000e+00
  [creep] predictor rel change = 0.000e+00

Convergence check


#### Iteration 2/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)


## Step 02/11: t = 1.00e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 1 | dt = 1.00e+07 s
Coupling = staggered
  → Max iterations              : 30
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.761e-01
  [creep] predictor rel change = 1.000e+00

Convergence check


#### Iteration 2/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.719e-01
  [creep] predictor rel change = 5.719e-01

Convergence check


#### Iteration 3/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 6.415e-02
  [creep] predictor rel change = 2.244e-01

Convergence check


#### Iteration 4/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.223e-02
  [creep] predictor rel change = 7.924e-02

Convergence check


#### Iteration 5/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 4.979e-04
  [creep] predictor rel change = 1.500e-02

Convergence check


#### Iteration 6/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 8.400e-07
  [creep] predictor rel change = 6.109e-04

Convergence check


#### Iteration 7/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.395e-12
  [creep] predictor rel change = 1.031e-06

Convergence check


#### Iteration 8/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 2.940e-12

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
  [creep] clad: max equivalent creep strain = 1.0081e-03
Computing symbolic result fields (strain, stress, ...)


## Step 03/11: t = 2.00e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 2 | dt = 1.00e+07 s
Coupling = staggered
  → Max iterations              : 30
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.315e-01
  [creep] predictor rel change = 2.323e+00

Convergence check


#### Iteration 2/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.058e-01
  [creep] predictor rel change = 5.719e-01

Convergence check


#### Iteration 3/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.948e-02
  [creep] predictor rel change = 2.244e-01

Convergence check


#### Iteration 4/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 7.526e-03
  [creep] predictor rel change = 7.924e-02

Convergence check


#### Iteration 5/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.064e-04
  [creep] predictor rel change = 1.500e-02

Convergence check


#### Iteration 6/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.170e-07
  [creep] predictor rel change = 6.109e-04

Convergence check


#### Iteration 7/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.473e-12
  [creep] predictor rel change = 1.031e-06

Convergence check


#### Iteration 8/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 2.938e-12

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
  [creep] clad: max equivalent creep strain = 2.0163e-03
Computing symbolic result fields (strain, stress, ...)


## Step 04/11: t = 3.00e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 3 | dt = 1.00e+07 s
Coupling = staggered
  → Max iterations              : 30
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.672e-01
  [creep] predictor rel change = 2.323e+00

Convergence check


#### Iteration 2/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 7.642e-02
  [creep] predictor rel change = 5.719e-01

Convergence check


#### Iteration 3/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.851e-02
  [creep] predictor rel change = 2.244e-01

Convergence check


#### Iteration 4/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.436e-03
  [creep] predictor rel change = 7.924e-02

Convergence check


#### Iteration 5/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.213e-04
  [creep] predictor rel change = 1.500e-02

Convergence check


#### Iteration 6/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.734e-07
  [creep] predictor rel change = 6.109e-04

Convergence check


#### Iteration 7/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.063e-12
  [creep] predictor rel change = 1.031e-06

Convergence check


#### Iteration 8/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 2.938e-12

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
  [creep] clad: max equivalent creep strain = 3.0244e-03
Computing symbolic result fields (strain, stress, ...)


## Step 05/11: t = 4.00e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 4 | dt = 1.00e+07 s
Coupling = staggered
  → Max iterations              : 30
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.308e-01
  [creep] predictor rel change = 2.323e+00

Convergence check


#### Iteration 2/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.981e-02
  [creep] predictor rel change = 5.719e-01

Convergence check


#### Iteration 3/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.231e-02
  [creep] predictor rel change = 2.244e-01

Convergence check


#### Iteration 4/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 4.254e-03
  [creep] predictor rel change = 7.924e-02

Convergence check


#### Iteration 5/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.732e-04
  [creep] predictor rel change = 1.500e-02

Convergence check


#### Iteration 6/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.922e-07
  [creep] predictor rel change = 6.109e-04

Convergence check


#### Iteration 7/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 8.322e-13
  [creep] predictor rel change = 1.031e-06

Convergence check


#### Iteration 8/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 2.937e-12

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
  [creep] clad: max equivalent creep strain = 4.0326e-03
Computing symbolic result fields (strain, stress, ...)


## Step 06/11: t = 5.00e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 5 | dt = 1.00e+07 s
Coupling = staggered
  → Max iterations              : 30
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.075e-01
  [creep] predictor rel change = 2.323e+00

Convergence check


#### Iteration 2/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 4.913e-02
  [creep] predictor rel change = 5.719e-01

Convergence check


#### Iteration 3/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.833e-02
  [creep] predictor rel change = 2.244e-01

Convergence check


#### Iteration 4/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.494e-03
  [creep] predictor rel change = 7.924e-02

Convergence check


#### Iteration 5/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.423e-04
  [creep] predictor rel change = 1.500e-02

Convergence check


#### Iteration 6/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.400e-07
  [creep] predictor rel change = 6.109e-04

Convergence check


#### Iteration 7/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 6.836e-13
  [creep] predictor rel change = 1.031e-06

Convergence check


#### Iteration 8/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 2.939e-12

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
  [creep] clad: max equivalent creep strain = 5.0407e-03
Computing symbolic result fields (strain, stress, ...)


## Step 07/11: t = 6.00e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 6 | dt = 1.00e+07 s
Coupling = staggered
  → Max iterations              : 30
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 9.119e-02
  [creep] predictor rel change = 2.323e+00

Convergence check


#### Iteration 2/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 4.169e-02
  [creep] predictor rel change = 5.719e-01

Convergence check


#### Iteration 3/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.555e-02
  [creep] predictor rel change = 2.244e-01

Convergence check


#### Iteration 4/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.965e-03
  [creep] predictor rel change = 7.924e-02

Convergence check


#### Iteration 5/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.207e-04
  [creep] predictor rel change = 1.500e-02

Convergence check


#### Iteration 6/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.037e-07
  [creep] predictor rel change = 6.109e-04

Convergence check


#### Iteration 7/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.789e-13
  [creep] predictor rel change = 1.031e-06

Convergence check


#### Iteration 8/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 2.938e-12

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
  [creep] clad: max equivalent creep strain = 6.0489e-03
Computing symbolic result fields (strain, stress, ...)


## Step 08/11: t = 7.00e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 7 | dt = 1.00e+07 s
Coupling = staggered
  → Max iterations              : 30
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 7.919e-02
  [creep] predictor rel change = 2.323e+00

Convergence check


#### Iteration 2/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.620e-02
  [creep] predictor rel change = 5.719e-01

Convergence check


#### Iteration 3/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.351e-02
  [creep] predictor rel change = 2.244e-01

Convergence check


#### Iteration 4/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.575e-03
  [creep] predictor rel change = 7.924e-02

Convergence check


#### Iteration 5/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.048e-04
  [creep] predictor rel change = 1.500e-02

Convergence check


#### Iteration 6/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.769e-07
  [creep] predictor rel change = 6.109e-04

Convergence check


#### Iteration 7/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.034e-13
  [creep] predictor rel change = 1.031e-06

Convergence check


#### Iteration 8/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 2.941e-12

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
  [creep] clad: max equivalent creep strain = 7.0570e-03
Computing symbolic result fields (strain, stress, ...)


## Step 09/11: t = 8.00e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 8 | dt = 1.00e+07 s
Coupling = staggered
  → Max iterations              : 30
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 6.998e-02
  [creep] predictor rel change = 2.323e+00

Convergence check


#### Iteration 2/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.199e-02
  [creep] predictor rel change = 5.719e-01

Convergence check


#### Iteration 3/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.193e-02
  [creep] predictor rel change = 2.244e-01

Convergence check


#### Iteration 4/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.275e-03
  [creep] predictor rel change = 7.924e-02

Convergence check


#### Iteration 5/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 9.264e-05
  [creep] predictor rel change = 1.500e-02

Convergence check


#### Iteration 6/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.563e-07
  [creep] predictor rel change = 6.109e-04

Convergence check


#### Iteration 7/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 4.453e-13
  [creep] predictor rel change = 1.031e-06

Convergence check


#### Iteration 8/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 2.944e-12

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
  [creep] clad: max equivalent creep strain = 8.0652e-03
Computing symbolic result fields (strain, stress, ...)


## Step 10/11: t = 9.00e+07 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 9 | dt = 1.00e+07 s
Coupling = staggered
  → Max iterations              : 30
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 6.269e-02
  [creep] predictor rel change = 2.323e+00

Convergence check


#### Iteration 2/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.866e-02
  [creep] predictor rel change = 5.719e-01

Convergence check


#### Iteration 3/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.069e-02
  [creep] predictor rel change = 2.244e-01

Convergence check


#### Iteration 4/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.038e-03
  [creep] predictor rel change = 7.924e-02

Convergence check


#### Iteration 5/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 8.299e-05
  [creep] predictor rel change = 1.500e-02

Convergence check


#### Iteration 6/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.400e-07
  [creep] predictor rel change = 6.109e-04

Convergence check


#### Iteration 7/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.986e-13
  [creep] predictor rel change = 1.031e-06

Convergence check


#### Iteration 8/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 2.942e-12

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
  [creep] clad: max equivalent creep strain = 9.0733e-03
Computing symbolic result fields (strain, stress, ...)


## Step 11/11: t = 1.00e+08 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 10 | dt = 1.00e+07 s
Coupling = staggered
  → Max iterations              : 30
  → Staggering tolerance |ΔT|   : 1.0e-08
  → Staggering tolerance |Δu|   : 1.0e-07
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-10
  → Relative tolerance mech     : 1.0e-09
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.678e-02
  [creep] predictor rel change = 2.323e+00

Convergence check


#### Iteration 2/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.596e-02
  [creep] predictor rel change = 5.719e-01

Convergence check


#### Iteration 3/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 9.683e-03
  [creep] predictor rel change = 2.244e-01

Convergence check


#### Iteration 4/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.846e-03
  [creep] predictor rel change = 7.924e-02

Convergence check


#### Iteration 5/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 7.516e-05
  [creep] predictor rel change = 1.500e-02

Convergence check


#### Iteration 6/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.268e-07
  [creep] predictor rel change = 6.109e-04

Convergence check


#### Iteration 7/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.611e-13
  [creep] predictor rel change = 1.031e-06

Convergence check


#### Iteration 8/30


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for clad, tag = 10
  → q_third[clad](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=600.00 K, max=600.00 K, mean=600.00 K
  ||ΔT||/||T|| = 0.000e+00

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating traction on region 4 → 100000000.0 Pa
  Building weak form, volume integrals (dx) for clad, tag = 10
  Applying mechanical traction on subdomain id = 4
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00
  [creep] predictor rel change = 2.948e-12

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
  [creep] clad: max equivalent creep strain = 1.0081e-02
Computing symbolic result fields (strain, stress, ...)

Simulation completed in 263.74 s
Total time steps solved: 11
