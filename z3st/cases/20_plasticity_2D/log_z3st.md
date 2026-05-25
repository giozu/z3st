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
  → Time steps          : 21
  → Regime              : 2d
  → Models active       :
      thermal    → OFF
      mechanical → ON
      damage     → OFF
      cluster    → OFF
      plasticity → ON
  → Gap conductance     : None (value = 0.0)



### FiniteElementSetup initializer

Mechanical element order: 1
Thermal function space (V_t): FunctionSpace(<Mesh #0>, Basix element (P, quadrilateral, 1, gll_warped, unset, False, float64, []))
Mechanical function space (V_m): FunctionSpace(<Mesh #0>, blocked element (Basix element (P, quadrilateral, 1, gll_warped, unset, False, float64, []), (2,)))
Scalar function space (Q): FunctionSpace(<Mesh #0>, Basix element (P, quadrilateral, 0, gll_warped, unset, True, float64, []))
Plasticity function spaces (V_pl_tensor, Q_pl) initializing with Quadrature degree 3
Plasticity function space (V_pl_tensor): FunctionSpace(<Mesh #0>, blocked element (QuadratureElement(quadrilateral, array([[0.21132487, 0.21132487],       [0.21132487, 0.78867513],       [0.78867513, 0.21132487],       [0.78867513, 0.78867513]]), array([0.25, 0.25, 0.25, 0.25]), IdentityPullback()), (3, 3)))
Plasticity function space (Q_pl): FunctionSpace(<Mesh #0>, QuadratureElement(quadrilateral, array([[0.21132487, 0.21132487],       [0.21132487, 0.78867513],       [0.78867513, 0.21132487],       [0.78867513, 0.78867513]]), array([0.25, 0.25, 0.25, 0.25]), IdentityPullback()))
[Solver] initializer
  Applied relaxation factor:
  → Temperature  : 0.9
  → Displacement : 1.0
  → Damage       : 0.4
  Adaptive relaxation disabled


[MechanicalModel] initializer
[MechanicalModel] options loaded from input.yaml:
  solver              : newton
  linear_solver       : iterative_hypre
  rtol                : 1e-06
  stag_tol            : 1e-06
  convergence         : rel_norm
  debug               : False
[PlasticityModel] initializer
[spine.load_materials]
Material loaded: steel
  → k defined as constant: 50.0
  → Gc not defined for steel
  → constitutive model: lame
  → constitutive model promoted to: plasticity (yield_strength present)
  E               → 200000000000.0 (float)
  G               → 76923076923.07692 (float)
  T_ref           → 300.0 (float)
  alpha           → 1e-05 (float)
  bulk_modulus    → 166666666666.66666 (float)
  constitutive_mode → plasticity (str)
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
  **[INFO]** Clamp_x mechanical BC on 'steel' → 0.0 (first step) at region 'xmax'
  **[INFO]** Clamp_x mechanical BC on 'steel' → 0.0 (first step) at region 'xmin'
  **[INFO]** Clamp_y mechanical BC on 'steel' → 0.0 (first step) at region 'ymin'
Computing symbolic result fields (strain, stress, ...)
[OutputWriter] interpolation not available for field on blocked element (Basix element (P, quadrilateral, 0, gll_warped, unset, True, float64, []), (3, 3)) (ValueError); falling back to L2 projection.
[OutputWriter] interpolation not available for field on Basix element (P, quadrilateral, 0, gll_warped, unset, True, float64, []) (ValueError); falling back to L2 projection.

**[INFO]** Hot-reload of allow-listed input.yaml parameters is active. Edit input.yaml during the run; changes apply at the next step boundary. Allowed keys: damage.{stag_tol,rtol,hybrid_constraint,gamma_star}, mechanical.{stag_tol,rtol}, thermal.{stag_tol,rtol}, solver_settings.{max_iters,relax_*}.


## Step 01/21: t = 0.00e+00 s | LHR = 0.00e+00 W/m

  → dt=0: solving static step / initial condition
Computing symbolic result fields (strain, stress, ...)



***


### spine - solve


***



Current step = 0 | dt = 0.00e+00 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 3 for integration measures.


#### Iteration 1/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 1 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)


## Step 02/21: t = 5.00e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 1 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 3 for integration measures.


#### Iteration 1/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 2e-05
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.000e+00

Convergence check


#### Iteration 2/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 2e-05
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.738e-09

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)


## Step 03/21: t = 1.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 2 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 3 for integration measures.


#### Iteration 1/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 4e-05
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.000e-01

Convergence check


#### Iteration 2/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 4e-05
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 8.439e-08

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)


## Step 04/21: t = 1.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 3 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 3 for integration measures.


#### Iteration 1/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 6e-05
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.333e-01

Convergence check


#### Iteration 2/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 6e-05
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.497e-08

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)


## Step 05/21: t = 2.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 4 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 3 for integration measures.


#### Iteration 1/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 8e-05
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.500e-01

Convergence check


#### Iteration 2/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 8e-05
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.627e-09

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)


## Step 06/21: t = 2.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 5 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 3 for integration measures.


#### Iteration 1/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0001
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.000e-01

Convergence check


#### Iteration 2/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0001
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.494e-08

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)


## Step 07/21: t = 3.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 6 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 3 for integration measures.


#### Iteration 1/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.00012
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.722e-01

Convergence check


#### Iteration 2/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.00012
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.673e-09

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)


## Step 08/21: t = 3.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 7 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 3 for integration measures.


#### Iteration 1/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.00014
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.281e-01

Convergence check


#### Iteration 2/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.00014
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.050e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)


## Step 09/21: t = 4.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 8 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 3 for integration measures.


#### Iteration 1/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.00016
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.883e-01

Convergence check


#### Iteration 2/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.00016
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 9.143e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)


## Step 10/21: t = 4.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 9 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 3 for integration measures.


#### Iteration 1/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.00018
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.602e-01

Convergence check


#### Iteration 2/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.00018
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 9.071e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)


## Step 11/21: t = 5.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 10 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 3 for integration measures.


#### Iteration 1/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0002
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.393e-01

Convergence check


#### Iteration 2/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0002
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.750e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)


## Step 12/21: t = 5.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 11 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 3 for integration measures.


#### Iteration 1/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.00022
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.232e-01

Convergence check


#### Iteration 2/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.00022
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.020e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)


## Step 13/21: t = 6.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 12 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 3 for integration measures.


#### Iteration 1/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.00024
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.103e-01

Convergence check


#### Iteration 2/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.00024
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 4.480e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)


## Step 14/21: t = 6.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 13 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 3 for integration measures.


#### Iteration 1/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.00026
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 9.982e-02

Convergence check


#### Iteration 2/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.00026
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.880e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)


## Step 15/21: t = 7.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 14 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 3 for integration measures.


#### Iteration 1/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.00028
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 9.111e-02

Convergence check


#### Iteration 2/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.00028
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.234e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)


## Step 16/21: t = 7.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 15 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 3 for integration measures.


#### Iteration 1/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0003
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 8.377e-02

Convergence check


#### Iteration 2/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0003
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.892e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)


## Step 17/21: t = 8.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 16 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 3 for integration measures.


#### Iteration 1/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.00032
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 7.749e-02

Convergence check


#### Iteration 2/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.00032
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.079e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)


## Step 18/21: t = 8.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 17 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 3 for integration measures.


#### Iteration 1/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.00034
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 7.207e-02

Convergence check


#### Iteration 2/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.00034
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.009e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)


## Step 19/21: t = 9.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 18 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 3 for integration measures.


#### Iteration 1/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.00036
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 6.734e-02

Convergence check


#### Iteration 2/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.00036
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.478e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)


## Step 20/21: t = 9.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 19 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 3 for integration measures.


#### Iteration 1/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.00038
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 6.318e-02

Convergence check


#### Iteration 2/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.00038
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.477e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)


## Step 21/21: t = 1.00e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 20 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 20
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 3 for integration measures.


#### Iteration 1/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0004
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.949e-02

Convergence check


#### Iteration 2/20


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0004
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  Building weak form, volume integrals (dx) for steel, tag = 5
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.286e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)

Simulation completed in 24.37 s
Total time steps solved: 21
