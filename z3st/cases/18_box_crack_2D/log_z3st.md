Info    : Reading 'mesh.msh'...
Info    : 13 entities
Info    : 42971 nodes
Info    : 86565 elements
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
  → Regime              : 2d
  → Models active       :
      thermal    → OFF
      mechanical → ON
      damage     → ON
      cluster    → OFF
      plasticity → OFF
      contact    → OFF
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
  → Displacement : 0.8
  → Damage       : 0.6
  Adaptive relaxation enabled
  → relax_growth  : 1.2
  → relax_shrink : 0.8
  → relax_min  : 0.05
  → relax_max : 0.95


[MechanicalModel] initializer
[MechanicalModel] options loaded from input.yaml:
  solver              : non-linear
  linear_solver       : iterative_hypre
  rtol                : 1e-06
  stag_tol            : 1e-06
  convergence         : rel_norm
DamageModel initializer
Options loaded from input.yaml:
  type                : AT2
  solver              : linear
  linear_solver       : iterative_hypre
  rtol                : 1e-06
  stag_tol            : 1e-06
  convergence         : rel_norm
  lc                  : 0.002
[spine.load_materials]
Material loaded: steel
  → k defined as constant: 45.0
  → Gc defined as constant: 2700.0
  - Material 'steel': sigma_c (AT2) from Gc = 2700.00 J/m2
  → constitutive model: lame
[spine.initialize_fields]

Initializing the displacement field...
  Initial u: min=0.00e+00 m, max=0.00e+00 m, mean=0.00e+00 m

Initializing the damage field...



***


### spine - set_boundary_conditions


***



Loading boundary conditions from 'boundary_conditions.yaml'
  **[INFO]** Clamp_x mechanical BC on 'steel' → 0.0 (first step) at region 'xmax'
  **[INFO]** Clamp_y mechanical BC on 'steel' → 0.0 (first step) at region 'ymin'
  **[INFO]** Neumann mechanical BC on 'steel' → ymax: 200000000.0 Pa (list loaded)

Setting damage boundary conditions...
  **[INFO]** Dirichlet damage BC on 'steel' → D = 1.0 at region 'crack'
Computing symbolic result fields (strain, stress, ...)

**[INFO]** Hot-reload of allow-listed input.yaml parameters is active. Edit input.yaml during the run; changes apply at the next step boundary. Allowed keys: damage.{stag_tol,rtol,hybrid_constraint,gamma_star}, mechanical.{stag_tol,rtol}, thermal.{stag_tol,rtol}, solver_settings.{max_iters,relax_*}.


## Step 01/1: t = 0.00e+00 s | LHR = 0.00e+00 W/m

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
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating traction on region 3 → 2.0e8 Pa
  Building weak form, volume integrals (dx) for steel, tag = 6
  Applying mechanical traction on subdomain id = 3
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.000e+00
  [adaptive] relax_u=0.80

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.000e+00
  [adaptive] relax_D=0.60
  |ΔD|_∞ = 6.000e-01

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating traction on region 3 → 2.0e8 Pa
  Building weak form, volume integrals (dx) for steel, tag = 6
  Applying mechanical traction on subdomain id = 3
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.000e-01
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.354e-01
  [adaptive] relax_D=0.72
  |ΔD|_∞ = 2.400e-01

Convergence check


#### Iteration 3/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating traction on region 3 → 2.0e8 Pa
  Building weak form, volume integrals (dx) for steel, tag = 6
  Applying mechanical traction on subdomain id = 3
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 4.750e-02
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.227e-01
  [adaptive] relax_D=0.86
  |ΔD|_∞ = 1.152e-01

Convergence check


#### Iteration 4/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating traction on region 3 → 2.0e8 Pa
  Building weak form, volume integrals (dx) for steel, tag = 6
  Applying mechanical traction on subdomain id = 3
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.375e-03
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.572e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.871e-02

Convergence check


#### Iteration 5/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating traction on region 3 → 2.0e8 Pa
  Building weak form, volume integrals (dx) for steel, tag = 6
  Applying mechanical traction on subdomain id = 3
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.188e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.137e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.788e-03

Convergence check


#### Iteration 6/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating traction on region 3 → 2.0e8 Pa
  Building weak form, volume integrals (dx) for steel, tag = 6
  Applying mechanical traction on subdomain id = 3
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.937e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.712e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.894e-04

Convergence check


#### Iteration 7/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating traction on region 3 → 2.0e8 Pa
  Building weak form, volume integrals (dx) for steel, tag = 6
  Applying mechanical traction on subdomain id = 3
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.969e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.869e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.447e-05

Convergence check


#### Iteration 8/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating traction on region 3 → 2.0e8 Pa
  Building weak form, volume integrals (dx) for steel, tag = 6
  Applying mechanical traction on subdomain id = 3
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.484e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.441e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.235e-07

Convergence check


#### Iteration 9/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 2 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 1 → 0.0
  **[INFO]** Updating traction on region 3 → 2.0e8 Pa
  Building weak form, volume integrals (dx) for steel, tag = 6
  Applying mechanical traction on subdomain id = 3
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 7.422e-10
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.236e-08
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.618e-08

Convergence check

**[SUCCESS]** Staggered solver converged in 9 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.4951e+04 J
  → Fracture energy : 4.0841e+03 J
  → Total energy    : 3.9035e+04 J

Simulation completed in 19.48 s
Total time steps solved: 1
