Info    : Reading 'mesh.msh'...
Info    : 75937 nodes
Info    : 152372 elements
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
  type                : AT1
  solver              : linear
  linear_solver       : iterative_hypre
  rtol                : 1e-05
  stag_tol            : 1e-05
  convergence         : rel_norm
  lc                  : 4e-06
[spine.load_materials]
Material loaded: uo2
  → k defined as constant: 5.0
  → Gc not defined for uo2
  - Material 'uo2': Gc (AT1) from sigma_c = 1.50e+08 Pa
  → constitutive model: lame
  E               → 358000000000.0 (float)
  G               → 145528455284.55286 (float)
  Gc              → 0.6703910614525139 (float)
  T_initial       → 1023.15 (float)
  T_ref           → 1023.15 (float)
  alpha           → 1e-05 (float)
  bulk_modulus    → 220987654320.98764 (float)
  constitutive_mode → lame (str)
  cp              → 280.0 (float)
  k               → 5.0 (float)
  lmbda           → 123968684131.28575 (float)
  name            → UO2 (str)
  nu              → 0.23 (float)
  rho             → 10970.0 (float)
  sigma_c         → 150000000.0 (float)
[spine.initialize_fields]

Initializing the displacement field...
  Initial u: min=0.00e+00 m, max=0.00e+00 m, mean=0.00e+00 m

Initializing the damage field...



***


### spine - set_boundary_conditions


***



Loading boundary conditions from 'boundary_conditions.yaml'
  **[INFO]** Constant Dirichlet vector (2D) → [0.0, 0.0]
  **[INFO]** Dirichlet mechanical BC on 'uo2' → [0.0, 0.0] at region 'ymin'
  **[INFO]** Step-dependent Dirichlet list (2D), length 21
  **[INFO]** Dirichlet mechanical BC on 'uo2' → [0.0, 0.0] at region 'ymax'

Setting damage boundary conditions...
  **[INFO]** Dirichlet damage BC on 'uo2' → D = 1.0 at region 'crack'
Computing symbolic result fields (strain, stress, ...)

**[INFO]** Hot-reload of allow-listed input.yaml parameters is active. Edit input.yaml during the run; changes apply at the next step boundary. Allowed keys: damage.{stag_tol,rtol,hybrid_constraint,gamma_star}, mechanical.{stag_tol,rtol}, thermal.{stag_tol,rtol}, solver_settings.{max_iters,relax_*}.


## Step 01/21: t = 0.00e+00 s | LHR = 0.00e+00 W/m

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
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 0.0]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.000e+00
  |ΔD|_∞ = 8.000e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 0.0]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 2.000e-01
  |ΔD|_∞ = 1.600e-01

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 0.0]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 4.000e-02
  |ΔD|_∞ = 3.200e-02

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 0.0]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 8.000e-03
  |ΔD|_∞ = 6.400e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 0.0]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.600e-03
  |ΔD|_∞ = 1.280e-03

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 0.0]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 3.200e-04
  |ΔD|_∞ = 2.560e-04

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 0.0]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 6.400e-05
  |ΔD|_∞ = 5.120e-05

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 0.0]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.280e-05
  |ΔD|_∞ = 1.024e-05

Convergence check


#### Iteration 9/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 0.0]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 2.560e-06
  |ΔD|_∞ = 2.048e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 9 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 0.0000e+00 J
  → Fracture energy : 1.3008e-03 J
  → Total energy    : 1.3008e-03 J


## Step 02/21: t = 1.80e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 1 | dt = 1.80e+02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5e-08]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 4.919e-02
  |ΔD|_∞ = 1.898e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5e-08]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.253e-18

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 9.837e-03
  |ΔD|_∞ = 3.797e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5e-08]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.253e-18

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.967e-03
  |ΔD|_∞ = 7.594e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5e-08]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 3.935e-04
  |ΔD|_∞ = 1.519e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5e-08]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.072e-16

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 7.870e-05
  |ΔD|_∞ = 3.037e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5e-08]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.074e-16

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.574e-05
  |ΔD|_∞ = 6.075e-05

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5e-08]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.747e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 3.148e-06
  |ΔD|_∞ = 1.215e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 7 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.8379e-04 J
  → Fracture energy : 1.2603e-03 J
  → Total energy    : 1.7441e-03 J


## Step 03/21: t = 3.60e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 2 | dt = 1.80e+02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.000e-01

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 3.946e-01
  |ΔD|_∞ = 4.449e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.904e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 7.892e-02
  |ΔD|_∞ = 8.898e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.324e-19

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.578e-02
  |ΔD|_∞ = 1.780e-02

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.135e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 3.157e-03
  |ΔD|_∞ = 3.559e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.136e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 6.314e-04
  |ΔD|_∞ = 7.118e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.963e-19

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.263e-04
  |ΔD|_∞ = 1.424e-04

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 2.525e-05
  |ΔD|_∞ = 2.847e-05

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.362e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 5.051e-06
  |ΔD|_∞ = 5.694e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9319e-03 J
  → Fracture energy : 7.1560e-04 J
  → Total energy    : 2.6475e-03 J


## Step 04/21: t = 5.40e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 3 | dt = 1.80e+02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.336e-01

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 2.764e-01
  |ΔD|_∞ = 3.379e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.446e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 5.528e-02
  |ΔD|_∞ = 6.759e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.446e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.106e-02
  |ΔD|_∞ = 1.352e-02

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.149e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 2.211e-03
  |ΔD|_∞ = 2.704e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.945e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 4.423e-04
  |ΔD|_∞ = 5.407e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.038e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 8.845e-05
  |ΔD|_∞ = 1.081e-04

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.769e-05
  |ΔD|_∞ = 2.163e-05

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 3.538e-06
  |ΔD|_∞ = 4.326e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.3176e-03 J
  → Fracture energy : 5.1544e-04 J
  → Total energy    : 4.8330e-03 J


## Step 05/21: t = 7.20e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 4 | dt = 1.80e+02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.509e-01

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 2.719e-01
  |ΔD|_∞ = 3.132e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.344e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 5.438e-02
  |ΔD|_∞ = 6.265e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.938e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.088e-02
  |ΔD|_∞ = 1.253e-02

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.952e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 2.175e-03
  |ΔD|_∞ = 2.506e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.952e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 4.350e-04
  |ΔD|_∞ = 5.012e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.932e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 8.701e-05
  |ΔD|_∞ = 1.002e-04

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.330e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.740e-05
  |ΔD|_∞ = 2.005e-05

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 3.480e-06
  |ΔD|_∞ = 4.009e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.5422e-03 J
  → Fracture energy : 4.1559e-04 J
  → Total energy    : 7.9578e-03 J


## Step 06/21: t = 9.00e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 5 | dt = 1.80e+02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.079e-01

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 2.599e-01
  |ΔD|_∞ = 4.464e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.288e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 5.198e-02
  |ΔD|_∞ = 8.927e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.288e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.040e-02
  |ΔD|_∞ = 1.785e-02

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.288e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 2.079e-03
  |ΔD|_∞ = 3.571e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 4.158e-04
  |ΔD|_∞ = 7.142e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 8.317e-05
  |ΔD|_∞ = 1.428e-04

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.663e-05
  |ΔD|_∞ = 2.857e-05

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.288e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 3.327e-06
  |ΔD|_∞ = 5.714e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0978e-02 J
  → Fracture energy : 4.1720e-04 J
  → Total energy    : 1.1395e-02 J


## Step 07/21: t = 1.08e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 6 | dt = 1.80e+02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.231e-01

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 3.518e-01
  |ΔD|_∞ = 7.216e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 7.035e-02
  |ΔD|_∞ = 1.443e-01

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.475e-16

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.407e-02
  |ΔD|_∞ = 2.887e-02

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.263e-16

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 2.814e-03
  |ΔD|_∞ = 5.773e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.262e-16

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 5.628e-04
  |ΔD|_∞ = 1.155e-03

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.113e-16

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.126e-04
  |ΔD|_∞ = 2.309e-04

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.041e-16

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 2.251e-05
  |ΔD|_∞ = 4.619e-05

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.475e-16

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 4.503e-06
  |ΔD|_∞ = 9.237e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0186e-02 J
  → Fracture energy : 5.6630e-04 J
  → Total energy    : 1.0752e-02 J


## Step 08/21: t = 1.26e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 7 | dt = 1.80e+02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.204e-01

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 3.947e-01
  |ΔD|_∞ = 5.287e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.622e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 7.895e-02
  |ΔD|_∞ = 1.057e-01

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.603e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.579e-02
  |ΔD|_∞ = 2.115e-02

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.610e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 3.158e-03
  |ΔD|_∞ = 4.229e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.036e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 6.316e-04
  |ΔD|_∞ = 8.459e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.263e-04
  |ΔD|_∞ = 1.692e-04

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 2.526e-05
  |ΔD|_∞ = 3.384e-05

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.036e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 5.053e-06
  |ΔD|_∞ = 6.767e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2166e-02 J
  → Fracture energy : 9.7704e-04 J
  → Total energy    : 1.3143e-02 J


## Step 09/21: t = 1.44e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 8 | dt = 1.80e+02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.468e-01

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 3.873e-01
  |ΔD|_∞ = 4.787e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.370e-16

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 7.746e-02
  |ΔD|_∞ = 9.574e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.318e-18

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.549e-02
  |ΔD|_∞ = 1.915e-02

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.370e-16

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 3.098e-03
  |ΔD|_∞ = 3.829e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.259e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 6.197e-04
  |ΔD|_∞ = 7.659e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.259e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.239e-04
  |ΔD|_∞ = 1.532e-04

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 2.479e-05
  |ΔD|_∞ = 3.064e-05

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 4.957e-06
  |ΔD|_∞ = 6.127e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3112e-02 J
  → Fracture energy : 2.0988e-03 J
  → Total energy    : 1.5211e-02 J


## Step 10/21: t = 1.62e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 9 | dt = 1.80e+02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.377e-01

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 3.396e-01
  |ΔD|_∞ = 4.145e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.663e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 6.792e-02
  |ΔD|_∞ = 8.290e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.649e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.358e-02
  |ΔD|_∞ = 1.658e-02

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 2.717e-03
  |ΔD|_∞ = 3.316e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.649e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 5.433e-04
  |ΔD|_∞ = 6.632e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.649e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.087e-04
  |ΔD|_∞ = 1.326e-04

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.649e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 2.173e-05
  |ΔD|_∞ = 2.653e-05

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.649e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 4.347e-06
  |ΔD|_∞ = 5.305e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0562e-02 J
  → Fracture energy : 4.8407e-03 J
  → Total energy    : 1.5403e-02 J


## Step 11/21: t = 1.80e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 10 | dt = 1.80e+02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.319e-01

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.899e-01
  |ΔD|_∞ = 3.423e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.148e-18

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 3.798e-02
  |ΔD|_∞ = 6.847e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.148e-18

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 7.596e-03
  |ΔD|_∞ = 1.369e-02

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.290e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.519e-03
  |ΔD|_∞ = 2.739e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 3.038e-04
  |ΔD|_∞ = 5.478e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.276e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 6.077e-05
  |ΔD|_∞ = 1.096e-04

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.148e-18

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.215e-05
  |ΔD|_∞ = 2.191e-05

Convergence check


#### Iteration 8/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.148e-18

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 2.431e-06
  |ΔD|_∞ = 4.382e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 8 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.2704e-03 J
  → Fracture energy : 7.3069e-03 J
  → Total energy    : 1.3577e-02 J


## Step 12/21: t = 1.98e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 11 | dt = 1.80e+02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.108e-01

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 5.658e-02
  |ΔD|_∞ = 3.045e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.132e-02
  |ΔD|_∞ = 6.090e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 2.263e-03
  |ΔD|_∞ = 1.218e-02

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.085e-16

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 4.526e-04
  |ΔD|_∞ = 2.436e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.091e-16

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 9.052e-05
  |ΔD|_∞ = 4.872e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.018e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.810e-05
  |ΔD|_∞ = 9.744e-05

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.092e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 3.621e-06
  |ΔD|_∞ = 1.949e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 7 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9881e-03 J
  → Fracture energy : 8.3983e-03 J
  → Total energy    : 1.1386e-02 J


## Step 13/21: t = 2.16e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 12 | dt = 1.80e+02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.123e-01

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.470e-02
  |ΔD|_∞ = 1.914e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 2.939e-03
  |ΔD|_∞ = 3.828e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.663e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 5.878e-04
  |ΔD|_∞ = 7.657e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.148e-18

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.176e-04
  |ΔD|_∞ = 1.531e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.148e-18

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 2.351e-05
  |ΔD|_∞ = 3.063e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.663e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 4.703e-06
  |ΔD|_∞ = 6.125e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.4739e-04 J
  → Fracture energy : 8.6133e-03 J
  → Total energy    : 8.9607e-03 J


## Step 14/21: t = 2.34e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 13 | dt = 1.80e+02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.255e-02

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.362e-03
  |ΔD|_∞ = 2.256e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.091e-18

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 2.724e-04
  |ΔD|_∞ = 4.513e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.091e-18

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 5.448e-05
  |ΔD|_∞ = 9.025e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.091e-18

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.090e-05
  |ΔD|_∞ = 1.805e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 2.179e-06
  |ΔD|_∞ = 3.610e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.8527e-05 J
  → Fracture energy : 8.6237e-03 J
  → Total energy    : 8.6823e-03 J


## Step 15/21: t = 2.52e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 14 | dt = 1.80e+02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.683e-02

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 3.915e-04
  |ΔD|_∞ = 1.161e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.374e-16

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 7.831e-05
  |ΔD|_∞ = 2.322e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.566e-05
  |ΔD|_∞ = 4.644e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.524e-16

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 3.132e-06
  |ΔD|_∞ = 9.289e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.3045e-05 J
  → Fracture energy : 8.6253e-03 J
  → Total energy    : 8.6883e-03 J


## Step 16/21: t = 2.70e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 15 | dt = 1.80e+02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.281e-02

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.428e-04
  |ΔD|_∞ = 7.111e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.498e-16

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 2.857e-05
  |ΔD|_∞ = 1.422e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.436e-16

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 5.714e-06
  |ΔD|_∞ = 2.845e-04

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.1316e-05 J
  → Fracture energy : 8.6255e-03 J
  → Total energy    : 8.6968e-03 J


## Step 17/21: t = 2.88e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 16 | dt = 1.80e+02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 8e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.948e-02

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 9.188e-05
  |ΔD|_∞ = 6.323e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 8e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.922e-18

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.838e-05
  |ΔD|_∞ = 1.265e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 8e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.822e-19

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 3.675e-06
  |ΔD|_∞ = 2.529e-04

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.0751e-05 J
  → Fracture energy : 8.6256e-03 J
  → Total energy    : 8.7063e-03 J


## Step 18/21: t = 3.06e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 17 | dt = 1.80e+02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 8.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.593e-02

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 6.365e-05
  |ΔD|_∞ = 5.284e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 8.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.233e-16

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.273e-05
  |ΔD|_∞ = 1.057e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 8.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.356e-16

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 2.546e-06
  |ΔD|_∞ = 2.114e-04

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.0768e-05 J
  → Fracture energy : 8.6257e-03 J
  → Total energy    : 8.7164e-03 J


## Step 19/21: t = 3.24e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 18 | dt = 1.80e+02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 9e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.571e-02

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 4.051e-05
  |ΔD|_∞ = 2.383e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 9e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 8.101e-06
  |ΔD|_∞ = 4.766e-04

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0141e-04 J
  → Fracture energy : 8.6257e-03 J
  → Total energy    : 8.7271e-03 J


## Step 20/21: t = 3.42e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 19 | dt = 1.80e+02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 9.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.515e-02

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 4.087e-05
  |ΔD|_∞ = 3.286e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 9.5e-07]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.367e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 8.174e-06
  |ΔD|_∞ = 6.573e-04

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1267e-04 J
  → Fracture energy : 8.6258e-03 J
  → Total energy    : 8.7384e-03 J


## Step 21/21: t = 3.60e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 20 | dt = 1.80e+02 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1e-06]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.329e-02

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 5.124e-05
  |ΔD|_∞ = 3.669e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1e-06]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.048e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 1.025e-05
  |ΔD|_∞ = 7.339e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1e-06]
  Building weak form, volume integrals (dx) for uo2, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.886e-17

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=6.70e-01, sigma_c=1.50e+08
  ||ΔD||/||D|| = 2.050e-06
  |ΔD|_∞ = 1.468e-04

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2444e-04 J
  → Fracture energy : 8.6258e-03 J
  → Total energy    : 8.7502e-03 J

Simulation completed in 286.99 s
Total time steps solved: 21
