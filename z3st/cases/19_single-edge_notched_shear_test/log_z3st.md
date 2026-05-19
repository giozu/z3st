Info    : Reading 'mesh.msh'...
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
  → Time steps          : 151
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
  linear_solver       : iterative_hypre
  rtol                : 1e-05
  stag_tol            : 1e-05
  convergence         : rel_norm
  lc                  : 4e-06
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
  **[INFO]** Step-dependent Dirichlet list (2D), length 151
  **[INFO]** Dirichlet mechanical BC on 'steel' → [0.0, 0.0] at region 'ymax'

Setting damage boundary conditions...
  **[INFO]** Dirichlet damage BC on 'steel' → D = 1.0 at region 'crack'
Computing symbolic result fields (strain, stress, ...)


## Step 01/151: t = 0.00e+00 s | LHR = 0.00e+00 W/m

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

**[SUCCESS]** Staggered solver converged in 9 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 0.0000e+00 J
  → Fracture energy : 1.3579e+00 J
  → Total energy    : 1.3579e+00 J


## Step 02/151: t = 2.40e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 1 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2e-07, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.131e-06
  |ΔD|_∞ = 3.008e-05

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2e-07, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.826e-06
  |ΔD|_∞ = 6.016e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0392e-03 J
  → Fracture energy : 1.3579e+00 J
  → Total energy    : 1.3589e+00 J


## Step 03/151: t = 4.80e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 2 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [4e-07, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.000e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.694e-05
  |ΔD|_∞ = 9.144e-05

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [4e-07, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.387e-06
  |ΔD|_∞ = 1.829e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.1566e-03 J
  → Fracture energy : 1.3579e+00 J
  → Total energy    : 1.3620e+00 J


## Step 04/151: t = 7.20e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 3 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6e-07, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.333e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.537e-05
  |ΔD|_∞ = 1.541e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6e-07, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.074e-06
  |ΔD|_∞ = 3.081e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.3522e-03 J
  → Fracture energy : 1.3579e+00 J
  → Total energy    : 1.3672e+00 J


## Step 05/151: t = 9.60e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 4 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8e-07, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.500e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.384e-05
  |ΔD|_∞ = 2.167e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8e-07, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.299e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.277e-05
  |ΔD|_∞ = 4.335e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8e-07, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.299e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.554e-06
  |ΔD|_∞ = 8.670e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6625e-02 J
  → Fracture energy : 1.3579e+00 J
  → Total energy    : 1.3745e+00 J


## Step 06/151: t = 1.20e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 5 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.000e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.031e-05
  |ΔD|_∞ = 2.725e-04

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
  ||ΔD||/||D|| = 1.606e-05
  |ΔD|_∞ = 5.451e-05

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
  ||ΔD||/||D|| = 3.212e-06
  |ΔD|_∞ = 1.090e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5976e-02 J
  → Fracture energy : 1.3579e+00 J
  → Total energy    : 1.3839e+00 J


## Step 07/151: t = 1.44e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 6 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.667e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.825e-05
  |ΔD|_∞ = 3.332e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.965e-05
  |ΔD|_∞ = 6.665e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.930e-06
  |ΔD|_∞ = 1.333e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.7403e-02 J
  → Fracture energy : 1.3579e+00 J
  → Total energy    : 1.3953e+00 J


## Step 08/151: t = 1.68e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 7 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.429e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.162e-04
  |ΔD|_∞ = 3.940e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.000e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.324e-05
  |ΔD|_∞ = 7.881e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.649e-06
  |ΔD|_∞ = 1.576e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.0906e-02 J
  → Fracture energy : 1.3579e+00 J
  → Total energy    : 1.4088e+00 J


## Step 09/151: t = 1.92e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 8 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.250e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.343e-04
  |ΔD|_∞ = 4.549e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.686e-05
  |ΔD|_∞ = 9.099e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.372e-06
  |ΔD|_∞ = 1.820e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.6483e-02 J
  → Fracture energy : 1.3579e+00 J
  → Total energy    : 1.4244e+00 J


## Step 10/151: t = 2.16e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 9 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.111e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.523e-04
  |ΔD|_∞ = 5.159e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.047e-05
  |ΔD|_∞ = 1.032e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.094e-06
  |ΔD|_∞ = 2.064e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.4134e-02 J
  → Fracture energy : 1.3579e+00 J
  → Total energy    : 1.4420e+00 J


## Step 11/151: t = 2.40e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 10 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.000e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.705e-04
  |ΔD|_∞ = 5.770e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.410e-05
  |ΔD|_∞ = 1.154e-04

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
  ||ΔD||/||D|| = 6.821e-06
  |ΔD|_∞ = 2.308e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0386e-01 J
  → Fracture energy : 1.3579e+00 J
  → Total energy    : 1.4618e+00 J


## Step 12/151: t = 2.64e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 11 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.092e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.889e-04
  |ΔD|_∞ = 6.382e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.373e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.778e-05
  |ΔD|_∞ = 1.276e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.555e-06
  |ΔD|_∞ = 2.553e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2565e-01 J
  → Fracture energy : 1.3579e+00 J
  → Total energy    : 1.4836e+00 J


## Step 13/151: t = 2.88e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 12 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.335e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.074e-04
  |ΔD|_∞ = 6.995e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.399e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.147e-05
  |ΔD|_∞ = 1.399e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.294e-06
  |ΔD|_∞ = 2.798e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4952e-01 J
  → Fracture energy : 1.3579e+00 J
  → Total energy    : 1.5074e+00 J


## Step 14/151: t = 3.12e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 13 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.694e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.258e-04
  |ΔD|_∞ = 7.610e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.740e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.516e-05
  |ΔD|_∞ = 1.522e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.602e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.032e-06
  |ΔD|_∞ = 3.044e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7545e-01 J
  → Fracture energy : 1.3579e+00 J
  → Total energy    : 1.5334e+00 J


## Step 15/151: t = 3.36e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 14 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.145e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.443e-04
  |ΔD|_∞ = 8.227e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.886e-05
  |ΔD|_∞ = 1.645e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.773e-06
  |ΔD|_∞ = 3.291e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0345e-01 J
  → Fracture energy : 1.3580e+00 J
  → Total energy    : 1.5614e+00 J


## Step 16/151: t = 3.60e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 15 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.669e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.632e-04
  |ΔD|_∞ = 8.846e-04

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
  ||ΔD||/||D|| = 5.264e-05
  |ΔD|_∞ = 1.769e-04

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
  ||ΔD||/||D|| = 1.053e-05
  |ΔD|_∞ = 3.538e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.105e-06
  |ΔD|_∞ = 7.077e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3351e-01 J
  → Fracture energy : 1.3580e+00 J
  → Total energy    : 1.5915e+00 J


## Step 17/151: t = 3.84e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 16 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.252e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.804e-04
  |ΔD|_∞ = 9.410e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.609e-05
  |ΔD|_∞ = 1.882e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.122e-05
  |ΔD|_∞ = 3.764e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.316e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.244e-06
  |ΔD|_∞ = 7.528e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6563e-01 J
  → Fracture energy : 1.3580e+00 J
  → Total energy    : 1.6236e+00 J


## Step 18/151: t = 4.08e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 17 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.885e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.995e-04
  |ΔD|_∞ = 1.003e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.990e-05
  |ΔD|_∞ = 2.006e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.198e-05
  |ΔD|_∞ = 4.011e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.396e-06
  |ΔD|_∞ = 8.023e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9981e-01 J
  → Fracture energy : 1.3580e+00 J
  → Total energy    : 1.6578e+00 J


## Step 19/151: t = 4.32e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 18 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.558e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.187e-04
  |ΔD|_∞ = 1.065e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.374e-05
  |ΔD|_∞ = 2.130e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.275e-05
  |ΔD|_∞ = 4.260e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.549e-06
  |ΔD|_∞ = 8.520e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3605e-01 J
  → Fracture energy : 1.3581e+00 J
  → Total energy    : 1.6941e+00 J


## Step 20/151: t = 4.56e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 19 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.266e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.381e-04
  |ΔD|_∞ = 1.127e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.763e-05
  |ΔD|_∞ = 2.255e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.992e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.353e-05
  |ΔD|_∞ = 4.509e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.992e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.705e-06
  |ΔD|_∞ = 9.019e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.7434e-01 J
  → Fracture energy : 1.3581e+00 J
  → Total energy    : 1.7324e+00 J


## Step 21/151: t = 4.80e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 20 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.003e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.577e-04
  |ΔD|_∞ = 1.190e-03

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
  ||ΔD||/||D|| = 7.155e-05
  |ΔD|_∞ = 2.380e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.671e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.431e-05
  |ΔD|_∞ = 4.760e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.671e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.862e-06
  |ΔD|_∞ = 9.520e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.1469e-01 J
  → Fracture energy : 1.3581e+00 J
  → Total energy    : 1.7728e+00 J


## Step 22/151: t = 5.04e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 21 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [4.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.765e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.780e-04
  |ΔD|_∞ = 1.253e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [4.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.561e-05
  |ΔD|_∞ = 2.506e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [4.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.512e-05
  |ΔD|_∞ = 5.012e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [4.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.024e-06
  |ΔD|_∞ = 1.002e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.5708e-01 J
  → Fracture energy : 1.3582e+00 J
  → Total energy    : 1.8153e+00 J


## Step 23/151: t = 5.28e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 22 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [4.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.549e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.980e-04
  |ΔD|_∞ = 1.316e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [4.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.959e-05
  |ΔD|_∞ = 2.632e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [4.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.592e-05
  |ΔD|_∞ = 5.265e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [4.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.184e-06
  |ΔD|_∞ = 1.053e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.0151e-01 J
  → Fracture energy : 1.3583e+00 J
  → Total energy    : 1.8598e+00 J


## Step 24/151: t = 5.52e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 23 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [4.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.352e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.190e-04
  |ΔD|_∞ = 1.380e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [4.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.308e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.381e-05
  |ΔD|_∞ = 2.759e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [4.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.676e-05
  |ΔD|_∞ = 5.519e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [4.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.132e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.352e-06
  |ΔD|_∞ = 1.104e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.4799e-01 J
  → Fracture energy : 1.3583e+00 J
  → Total energy    : 1.9063e+00 J


## Step 25/151: t = 5.76e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 24 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [4.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.171e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.392e-04
  |ΔD|_∞ = 1.444e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [4.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.783e-05
  |ΔD|_∞ = 2.887e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [4.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.757e-05
  |ΔD|_∞ = 5.774e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [4.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.513e-06
  |ΔD|_∞ = 1.155e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.9650e-01 J
  → Fracture energy : 1.3584e+00 J
  → Total energy    : 1.9549e+00 J


## Step 26/151: t = 6.00e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 25 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.004e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.619e-04
  |ΔD|_∞ = 1.508e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.238e-05
  |ΔD|_∞ = 3.016e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.848e-05
  |ΔD|_∞ = 6.032e-05

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
  ||ΔD||/||D|| = 3.695e-06
  |ΔD|_∞ = 1.206e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.4704e-01 J
  → Fracture energy : 1.3585e+00 J
  → Total energy    : 2.0056e+00 J


## Step 27/151: t = 6.24e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 26 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.851e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.888e-04
  |ΔD|_∞ = 1.799e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.775e-05
  |ΔD|_∞ = 3.597e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.955e-05
  |ΔD|_∞ = 7.195e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.910e-06
  |ΔD|_∞ = 1.439e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.9962e-01 J
  → Fracture energy : 1.3586e+00 J
  → Total energy    : 2.0582e+00 J


## Step 28/151: t = 6.48e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 27 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.708e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.039e-04
  |ΔD|_∞ = 1.638e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.008e-04
  |ΔD|_∞ = 3.275e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.016e-05
  |ΔD|_∞ = 6.550e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.031e-06
  |ΔD|_∞ = 1.310e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.5421e-01 J
  → Fracture energy : 1.3587e+00 J
  → Total energy    : 2.1130e+00 J


## Step 29/151: t = 6.72e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 28 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.576e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.373e-04
  |ΔD|_∞ = 2.704e-03

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
  ||ΔD||/||D|| = 1.075e-04
  |ΔD|_∞ = 5.408e-04

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
  ||ΔD||/||D|| = 2.149e-05
  |ΔD|_∞ = 1.082e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.298e-06
  |ΔD|_∞ = 2.163e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.1082e-01 J
  → Fracture energy : 1.3589e+00 J
  → Total energy    : 2.1697e+00 J


## Step 30/151: t = 6.96e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 29 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.454e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.627e-04
  |ΔD|_∞ = 2.925e-03

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
  ||ΔD||/||D|| = 1.125e-04
  |ΔD|_∞ = 5.851e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.251e-05
  |ΔD|_∞ = 1.170e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [5.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.200e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.501e-06
  |ΔD|_∞ = 2.340e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.6945e-01 J
  → Fracture energy : 1.3590e+00 J
  → Total energy    : 2.2285e+00 J


## Step 31/151: t = 7.20e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 30 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.339e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.768e-04
  |ΔD|_∞ = 1.835e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.154e-04
  |ΔD|_∞ = 3.671e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.307e-05
  |ΔD|_∞ = 7.341e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.614e-06
  |ΔD|_∞ = 1.468e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.3008e-01 J
  → Fracture energy : 1.3592e+00 J
  → Total energy    : 2.2893e+00 J


## Step 32/151: t = 7.44e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 31 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.232e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.989e-04
  |ΔD|_∞ = 1.902e-03

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
  ||ΔD||/||D|| = 1.198e-04
  |ΔD|_∞ = 3.805e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.223e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.396e-05
  |ΔD|_∞ = 7.609e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.223e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.791e-06
  |ΔD|_∞ = 1.522e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.9271e-01 J
  → Fracture energy : 1.3594e+00 J
  → Total energy    : 2.3521e+00 J


## Step 33/151: t = 7.68e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 32 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.131e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.259e-04
  |ΔD|_∞ = 2.215e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.252e-04
  |ΔD|_∞ = 4.430e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.504e-05
  |ΔD|_∞ = 8.861e-05

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
  ||ΔD||/||D|| = 5.007e-06
  |ΔD|_∞ = 1.772e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0573e+00 J
  → Fracture energy : 1.3596e+00 J
  → Total energy    : 2.4169e+00 J


## Step 34/151: t = 7.92e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 33 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.037e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.485e-04
  |ΔD|_∞ = 2.038e-03

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
  ||ΔD||/||D|| = 1.297e-04
  |ΔD|_∞ = 4.075e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.594e-05
  |ΔD|_∞ = 8.151e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.188e-06
  |ΔD|_∞ = 1.630e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1240e+00 J
  → Fracture energy : 1.3598e+00 J
  → Total energy    : 2.4838e+00 J


## Step 35/151: t = 8.16e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 34 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.948e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.762e-04
  |ΔD|_∞ = 2.106e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [6.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.370e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.352e-04
  |ΔD|_∞ = 4.212e-04

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
  ||ΔD||/||D|| = 2.705e-05
  |ΔD|_∞ = 8.425e-05

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
  ||ΔD||/||D|| = 5.409e-06
  |ΔD|_∞ = 1.685e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1926e+00 J
  → Fracture energy : 1.3601e+00 J
  → Total energy    : 2.5527e+00 J


## Step 36/151: t = 8.40e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 35 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.864e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.161e-04
  |ΔD|_∞ = 2.386e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.432e-04
  |ΔD|_∞ = 4.771e-04

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
  ||ΔD||/||D|| = 2.864e-05
  |ΔD|_∞ = 9.542e-05

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
  ||ΔD||/||D|| = 5.729e-06
  |ΔD|_∞ = 1.908e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2631e+00 J
  → Fracture energy : 1.3604e+00 J
  → Total energy    : 2.6235e+00 J


## Step 37/151: t = 8.64e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 36 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.785e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.359e-04
  |ΔD|_∞ = 2.245e-03

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
  ||ΔD||/||D|| = 1.472e-04
  |ΔD|_∞ = 4.490e-04

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
  ||ΔD||/||D|| = 2.944e-05
  |ΔD|_∞ = 8.979e-05

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
  ||ΔD||/||D|| = 5.887e-06
  |ΔD|_∞ = 1.796e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3357e+00 J
  → Fracture energy : 1.3607e+00 J
  → Total energy    : 2.6964e+00 J


## Step 38/151: t = 8.88e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 37 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.711e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.645e-04
  |ΔD|_∞ = 2.315e-03

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
  ||ΔD||/||D|| = 1.529e-04
  |ΔD|_∞ = 4.631e-04

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
  ||ΔD||/||D|| = 3.058e-05
  |ΔD|_∞ = 9.261e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.116e-06
  |ΔD|_∞ = 1.852e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4101e+00 J
  → Fracture energy : 1.3610e+00 J
  → Total energy    : 2.7712e+00 J


## Step 39/151: t = 9.12e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 38 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.640e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.981e-04
  |ΔD|_∞ = 3.350e-03

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
  ||ΔD||/||D|| = 1.596e-04
  |ΔD|_∞ = 6.700e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.193e-05
  |ΔD|_∞ = 1.340e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.102e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.385e-06
  |ΔD|_∞ = 2.680e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4866e+00 J
  → Fracture energy : 1.3614e+00 J
  → Total energy    : 2.8480e+00 J


## Step 40/151: t = 9.36e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 39 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.573e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.813e-04
  |ΔD|_∞ = 6.411e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [7.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.763e-04
  |ΔD|_∞ = 1.282e-03

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
  ||ΔD||/||D|| = 3.525e-05
  |ΔD|_∞ = 2.565e-04

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
  ||ΔD||/||D|| = 7.050e-06
  |ΔD|_∞ = 5.129e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5649e+00 J
  → Fracture energy : 1.3618e+00 J
  → Total energy    : 2.9267e+00 J


## Step 41/151: t = 9.60e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 40 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.510e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.074e-04
  |ΔD|_∞ = 4.112e-03

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
  ||ΔD||/||D|| = 1.815e-04
  |ΔD|_∞ = 8.225e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.629e-05
  |ΔD|_∞ = 1.645e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.259e-06
  |ΔD|_∞ = 3.290e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6452e+00 J
  → Fracture energy : 1.3623e+00 J
  → Total energy    : 3.0075e+00 J


## Step 42/151: t = 9.84e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 41 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.449e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.183e-04
  |ΔD|_∞ = 3.523e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.837e-04
  |ΔD|_∞ = 7.046e-04

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
  ||ΔD||/||D|| = 3.673e-05
  |ΔD|_∞ = 1.409e-04

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
  ||ΔD||/||D|| = 7.346e-06
  |ΔD|_∞ = 2.818e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7273e+00 J
  → Fracture energy : 1.3628e+00 J
  → Total energy    : 3.0901e+00 J


## Step 43/151: t = 1.01e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 42 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.391e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.964e-04
  |ΔD|_∞ = 7.236e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.749e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.993e-04
  |ΔD|_∞ = 1.447e-03

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
  ||ΔD||/||D|| = 3.985e-05
  |ΔD|_∞ = 2.895e-04

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
  ||ΔD||/||D|| = 7.971e-06
  |ΔD|_∞ = 5.789e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8113e+00 J
  → Fracture energy : 1.3634e+00 J
  → Total energy    : 3.1747e+00 J


## Step 44/151: t = 1.03e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 43 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.337e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.041e-03
  |ΔD|_∞ = 5.557e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.610e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.082e-04
  |ΔD|_∞ = 1.111e-03

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
  ||ΔD||/||D|| = 4.164e-05
  |ΔD|_∞ = 2.223e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.911e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.327e-06
  |ΔD|_∞ = 4.445e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8972e+00 J
  → Fracture energy : 1.3640e+00 J
  → Total energy    : 3.2612e+00 J


## Step 45/151: t = 1.06e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 44 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.285e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.131e-03
  |ΔD|_∞ = 1.022e-02

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
  ||ΔD||/||D|| = 2.262e-04
  |ΔD|_∞ = 2.045e-03

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
  ||ΔD||/||D|| = 4.525e-05
  |ΔD|_∞ = 4.089e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [8.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.050e-06
  |ΔD|_∞ = 8.179e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9849e+00 J
  → Fracture energy : 1.3647e+00 J
  → Total energy    : 3.3496e+00 J


## Step 46/151: t = 1.08e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 45 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.236e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.135e-03
  |ΔD|_∞ = 5.011e-03

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
  ||ΔD||/||D|| = 2.270e-04
  |ΔD|_∞ = 1.002e-03

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
  ||ΔD||/||D|| = 4.541e-05
  |ΔD|_∞ = 2.004e-04

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
  ||ΔD||/||D|| = 9.081e-06
  |ΔD|_∞ = 4.009e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0745e+00 J
  → Fracture energy : 1.3654e+00 J
  → Total energy    : 3.4399e+00 J


## Step 47/151: t = 1.10e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 46 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.188e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.219e-03
  |ΔD|_∞ = 9.156e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.2e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.438e-04
  |ΔD|_∞ = 1.831e-03

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
  ||ΔD||/||D|| = 4.877e-05
  |ΔD|_∞ = 3.662e-04

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
  ||ΔD||/||D|| = 9.753e-06
  |ΔD|_∞ = 7.324e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1658e+00 J
  → Fracture energy : 1.3662e+00 J
  → Total energy    : 3.5320e+00 J


## Step 48/151: t = 1.13e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 47 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.4e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.143e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.252e-03
  |ΔD|_∞ = 3.235e-03

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
  ||ΔD||/||D|| = 2.505e-04
  |ΔD|_∞ = 6.470e-04

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
  ||ΔD||/||D|| = 5.009e-05
  |ΔD|_∞ = 1.294e-04

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
  ||ΔD||/||D|| = 1.002e-05
  |ΔD|_∞ = 2.588e-05

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
  ||ΔD||/||D|| = 2.004e-06
  |ΔD|_∞ = 5.176e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2589e+00 J
  → Fracture energy : 1.3671e+00 J
  → Total energy    : 3.6260e+00 J


## Step 49/151: t = 1.15e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 48 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.100e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.335e-03
  |ΔD|_∞ = 5.373e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.169e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.671e-04
  |ΔD|_∞ = 1.075e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.169e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.341e-05
  |ΔD|_∞ = 2.149e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.068e-05
  |ΔD|_∞ = 4.299e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.6e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.137e-06
  |ΔD|_∞ = 8.597e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3537e+00 J
  → Fracture energy : 1.3681e+00 J
  → Total energy    : 3.7218e+00 J


## Step 50/151: t = 1.18e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 49 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.059e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.433e-03
  |ΔD|_∞ = 6.863e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.516e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.867e-04
  |ΔD|_∞ = 1.373e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.516e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.734e-05
  |ΔD|_∞ = 2.745e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.147e-05
  |ΔD|_∞ = 5.490e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [9.8e-06, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.294e-06
  |ΔD|_∞ = 1.098e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4502e+00 J
  → Fracture energy : 1.3692e+00 J
  → Total energy    : 3.8194e+00 J


## Step 51/151: t = 1.20e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 50 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.020e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.538e-03
  |ΔD|_∞ = 6.626e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.077e-04
  |ΔD|_∞ = 1.325e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.153e-05
  |ΔD|_∞ = 2.651e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.231e-05
  |ΔD|_∞ = 5.301e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.461e-06
  |ΔD|_∞ = 1.060e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5482e+00 J
  → Fracture energy : 1.3704e+00 J
  → Total energy    : 3.9186e+00 J


## Step 52/151: t = 1.22e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 51 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.02e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.984e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.654e-03
  |ΔD|_∞ = 7.152e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.02e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.309e-04
  |ΔD|_∞ = 1.430e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.02e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.617e-05
  |ΔD|_∞ = 2.861e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.02e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.323e-05
  |ΔD|_∞ = 5.722e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.02e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.647e-06
  |ΔD|_∞ = 1.144e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6477e+00 J
  → Fracture energy : 1.3718e+00 J
  → Total energy    : 4.0195e+00 J


## Step 53/151: t = 1.25e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 52 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.04e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.949e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.814e-03
  |ΔD|_∞ = 9.321e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.04e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.621e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.628e-04
  |ΔD|_∞ = 1.864e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.04e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.255e-05
  |ΔD|_∞ = 3.729e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.04e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.451e-05
  |ΔD|_∞ = 7.457e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.04e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.902e-06
  |ΔD|_∞ = 1.491e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7487e+00 J
  → Fracture energy : 1.3734e+00 J
  → Total energy    : 4.1220e+00 J


## Step 54/151: t = 1.27e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 53 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.06e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.917e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.071e-03
  |ΔD|_∞ = 1.392e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.06e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.142e-04
  |ΔD|_∞ = 2.783e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.06e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.283e-05
  |ΔD|_∞ = 5.566e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.06e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.657e-05
  |ΔD|_∞ = 1.113e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.06e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.313e-06
  |ΔD|_∞ = 2.227e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8507e+00 J
  → Fracture energy : 1.3752e+00 J
  → Total energy    : 4.2259e+00 J


## Step 55/151: t = 1.30e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 54 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.08e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.889e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.252e-03
  |ΔD|_∞ = 1.017e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.08e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.504e-04
  |ΔD|_∞ = 2.033e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.08e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.008e-05
  |ΔD|_∞ = 4.067e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.08e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.802e-05
  |ΔD|_∞ = 8.134e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.08e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.108e-20

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.603e-06
  |ΔD|_∞ = 1.627e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9537e+00 J
  → Fracture energy : 1.3773e+00 J
  → Total energy    : 4.3310e+00 J


## Step 56/151: t = 1.32e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 55 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.1e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.863e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.539e-03
  |ΔD|_∞ = 1.541e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.1e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.078e-04
  |ΔD|_∞ = 3.081e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.1e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.016e-04
  |ΔD|_∞ = 6.163e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.1e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.057e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.031e-05
  |ΔD|_∞ = 1.233e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.1e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.457e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.062e-06
  |ΔD|_∞ = 2.465e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.0574e+00 J
  → Fracture energy : 1.3797e+00 J
  → Total energy    : 4.4371e+00 J


## Step 57/151: t = 1.34e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 56 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.12e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.841e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.032e-03
  |ΔD|_∞ = 1.629e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.12e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.063e-04
  |ΔD|_∞ = 3.258e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.12e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.213e-04
  |ΔD|_∞ = 6.515e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.12e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.425e-05
  |ΔD|_∞ = 1.303e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.12e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.851e-06
  |ΔD|_∞ = 2.606e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1610e+00 J
  → Fracture energy : 1.3828e+00 J
  → Total energy    : 4.5438e+00 J


## Step 58/151: t = 1.37e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 57 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.14e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.831e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.396e-03
  |ΔD|_∞ = 1.681e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.14e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.793e-04
  |ΔD|_∞ = 3.363e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.14e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.359e-04
  |ΔD|_∞ = 6.726e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.14e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.717e-05
  |ΔD|_∞ = 1.345e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.14e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.434e-06
  |ΔD|_∞ = 2.690e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.2642e+00 J
  → Fracture energy : 1.3864e+00 J
  → Total energy    : 4.6507e+00 J


## Step 59/151: t = 1.39e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 58 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.16e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.826e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.993e-03
  |ΔD|_∞ = 2.117e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.16e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.511e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.985e-04
  |ΔD|_∞ = 4.233e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.16e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.511e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.597e-04
  |ΔD|_∞ = 8.466e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.16e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.194e-05
  |ΔD|_∞ = 1.693e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.16e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.388e-06
  |ΔD|_∞ = 3.386e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3662e+00 J
  → Fracture energy : 1.3910e+00 J
  → Total energy    : 4.7572e+00 J


## Step 60/151: t = 1.42e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 59 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.18e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.839e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.598e-03
  |ΔD|_∞ = 2.208e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.18e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.153e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.196e-04
  |ΔD|_∞ = 4.415e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.18e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.165e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.839e-04
  |ΔD|_∞ = 8.831e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.18e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.165e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.679e-05
  |ΔD|_∞ = 1.766e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.18e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.357e-06
  |ΔD|_∞ = 3.532e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.4659e+00 J
  → Fracture energy : 1.3965e+00 J
  → Total energy    : 4.8624e+00 J


## Step 61/151: t = 1.44e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 60 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.2e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.873e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.205e-03
  |ΔD|_∞ = 2.107e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.2e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.041e-03
  |ΔD|_∞ = 4.213e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.2e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.082e-04
  |ΔD|_∞ = 8.427e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.2e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.164e-05
  |ΔD|_∞ = 1.685e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.2e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.328e-06
  |ΔD|_∞ = 3.371e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.5624e+00 J
  → Fracture energy : 1.4032e+00 J
  → Total energy    : 4.9656e+00 J


## Step 62/151: t = 1.46e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 61 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.22e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.927e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.000e-03
  |ΔD|_∞ = 2.480e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.22e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.200e-03
  |ΔD|_∞ = 4.961e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.22e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.307e-20

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.400e-04
  |ΔD|_∞ = 9.922e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.22e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.307e-20

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.800e-05
  |ΔD|_∞ = 1.984e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.22e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.871e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.600e-06
  |ΔD|_∞ = 3.969e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.6543e+00 J
  → Fracture energy : 1.4113e+00 J
  → Total energy    : 5.0657e+00 J


## Step 63/151: t = 1.49e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 62 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.24e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.030e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.774e-03
  |ΔD|_∞ = 3.095e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.24e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.355e-03
  |ΔD|_∞ = 6.190e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.24e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.710e-04
  |ΔD|_∞ = 1.238e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.24e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.419e-05
  |ΔD|_∞ = 2.476e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.24e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.084e-05
  |ΔD|_∞ = 4.952e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.24e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.168e-06
  |ΔD|_∞ = 9.904e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.7410e+00 J
  → Fracture energy : 1.4209e+00 J
  → Total energy    : 5.1619e+00 J


## Step 64/151: t = 1.51e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 63 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.26e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.153e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.375e-03
  |ΔD|_∞ = 2.821e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.26e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.475e-03
  |ΔD|_∞ = 5.642e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.26e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.950e-04
  |ΔD|_∞ = 1.128e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.26e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.900e-05
  |ΔD|_∞ = 2.257e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.26e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.788e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.180e-05
  |ΔD|_∞ = 4.513e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.26e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.788e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.360e-06
  |ΔD|_∞ = 9.027e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.8230e+00 J
  → Fracture energy : 1.4316e+00 J
  → Total energy    : 5.2545e+00 J


## Step 65/151: t = 1.54e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 64 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.28e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.268e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.128e-03
  |ΔD|_∞ = 3.689e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.28e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.626e-03
  |ΔD|_∞ = 7.377e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.28e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.251e-04
  |ΔD|_∞ = 1.475e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.28e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.502e-05
  |ΔD|_∞ = 2.951e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.28e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.300e-05
  |ΔD|_∞ = 5.902e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.28e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.687e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.601e-06
  |ΔD|_∞ = 1.180e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.9013e+00 J
  → Fracture energy : 1.4433e+00 J
  → Total energy    : 5.3446e+00 J


## Step 66/151: t = 1.56e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 65 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.3e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.405e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.080e-03
  |ΔD|_∞ = 3.358e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.3e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.181e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.816e-03
  |ΔD|_∞ = 6.716e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.3e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.181e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.632e-04
  |ΔD|_∞ = 1.343e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.3e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.264e-05
  |ΔD|_∞ = 2.686e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.3e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.453e-05
  |ΔD|_∞ = 5.373e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.3e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.905e-06
  |ΔD|_∞ = 1.075e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.9748e+00 J
  → Fracture energy : 1.4562e+00 J
  → Total energy    : 5.4310e+00 J


## Step 67/151: t = 1.58e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 66 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.32e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.676e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.813e-03
  |ΔD|_∞ = 3.493e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.32e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.015e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.963e-03
  |ΔD|_∞ = 6.986e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.32e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.015e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.925e-04
  |ΔD|_∞ = 1.397e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.32e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.850e-05
  |ΔD|_∞ = 2.794e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.32e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.570e-05
  |ΔD|_∞ = 5.589e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.32e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.140e-06
  |ΔD|_∞ = 1.118e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.0411e+00 J
  → Fracture energy : 1.4703e+00 J
  → Total energy    : 5.5113e+00 J


## Step 68/151: t = 1.61e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 67 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.34e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.004e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.047e-02
  |ΔD|_∞ = 4.759e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.34e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.094e-03
  |ΔD|_∞ = 9.519e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.34e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.187e-04
  |ΔD|_∞ = 1.904e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.34e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.174e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.374e-05
  |ΔD|_∞ = 3.807e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.34e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.174e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.675e-05
  |ΔD|_∞ = 7.615e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.34e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.350e-06
  |ΔD|_∞ = 1.523e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.0979e+00 J
  → Fracture energy : 1.4852e+00 J
  → Total energy    : 5.5830e+00 J


## Step 69/151: t = 1.63e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 68 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.36e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.376e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.036e-02
  |ΔD|_∞ = 5.472e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.36e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.071e-03
  |ΔD|_∞ = 1.094e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.36e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.142e-04
  |ΔD|_∞ = 2.189e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.36e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.285e-05
  |ΔD|_∞ = 4.377e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.36e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.657e-05
  |ΔD|_∞ = 8.755e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.36e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.314e-06
  |ΔD|_∞ = 1.751e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.1595e+00 J
  → Fracture energy : 1.4987e+00 J
  → Total energy    : 5.6583e+00 J


## Step 70/151: t = 1.66e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 69 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.38e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.760e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.799e-03
  |ΔD|_∞ = 6.527e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.38e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.960e-03
  |ΔD|_∞ = 1.305e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.38e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.919e-04
  |ΔD|_∞ = 2.611e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.38e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.839e-05
  |ΔD|_∞ = 5.222e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.38e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.568e-05
  |ΔD|_∞ = 1.044e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.38e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.135e-06
  |ΔD|_∞ = 2.089e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.2447e+00 J
  → Fracture energy : 1.5093e+00 J
  → Total energy    : 5.7541e+00 J


## Step 71/151: t = 1.68e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 70 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.4e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.900e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.545e-03
  |ΔD|_∞ = 7.289e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.4e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.909e-03
  |ΔD|_∞ = 1.458e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.4e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.818e-04
  |ΔD|_∞ = 2.915e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.4e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.636e-05
  |ΔD|_∞ = 5.831e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.4e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.527e-05
  |ΔD|_∞ = 1.166e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.4e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.054e-06
  |ΔD|_∞ = 2.332e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.3487e+00 J
  → Fracture energy : 1.5181e+00 J
  → Total energy    : 5.8668e+00 J


## Step 72/151: t = 1.70e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 71 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.42e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.575e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.796e-03
  |ΔD|_∞ = 8.096e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.42e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.959e-03
  |ΔD|_∞ = 1.619e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.42e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.918e-04
  |ΔD|_∞ = 3.238e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.42e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.837e-05
  |ΔD|_∞ = 6.476e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.42e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.567e-05
  |ΔD|_∞ = 1.295e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.42e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.135e-06
  |ΔD|_∞ = 2.591e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.4569e+00 J
  → Fracture energy : 1.5269e+00 J
  → Total energy    : 5.9838e+00 J


## Step 73/151: t = 1.73e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 72 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.44e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.512e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.025e-02
  |ΔD|_∞ = 8.858e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.44e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.218e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.051e-03
  |ΔD|_∞ = 1.772e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.44e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.218e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.101e-04
  |ΔD|_∞ = 3.543e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.44e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.203e-05
  |ΔD|_∞ = 7.087e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.44e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.641e-05
  |ΔD|_∞ = 1.417e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.44e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.281e-06
  |ΔD|_∞ = 2.835e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.5659e+00 J
  → Fracture energy : 1.5366e+00 J
  → Total energy    : 6.1025e+00 J


## Step 74/151: t = 1.75e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 73 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.46e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.485e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.070e-02
  |ΔD|_∞ = 9.623e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.46e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.140e-03
  |ΔD|_∞ = 1.925e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.46e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.279e-04
  |ΔD|_∞ = 3.849e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.46e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.559e-05
  |ΔD|_∞ = 7.699e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.46e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.712e-05
  |ΔD|_∞ = 1.540e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.46e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.424e-06
  |ΔD|_∞ = 3.079e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.6736e+00 J
  → Fracture energy : 1.5472e+00 J
  → Total energy    : 6.2208e+00 J


## Step 75/151: t = 1.78e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 74 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.48e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.482e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.096e-02
  |ΔD|_∞ = 1.020e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.48e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.193e-03
  |ΔD|_∞ = 2.039e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.48e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.790e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.386e-04
  |ΔD|_∞ = 4.079e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.48e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.790e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.772e-05
  |ΔD|_∞ = 8.158e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.48e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.754e-05
  |ΔD|_∞ = 1.632e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.48e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.509e-06
  |ΔD|_∞ = 3.263e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.7802e+00 J
  → Fracture energy : 1.5588e+00 J
  → Total energy    : 6.3390e+00 J


## Step 76/151: t = 1.80e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 75 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.5e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.478e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.096e-02
  |ΔD|_∞ = 1.059e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.5e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.192e-03
  |ΔD|_∞ = 2.119e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.5e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.384e-04
  |ΔD|_∞ = 4.238e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.5e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.765e-05
  |ΔD|_∞ = 8.476e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.5e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.753e-05
  |ΔD|_∞ = 1.695e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.5e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.506e-06
  |ΔD|_∞ = 3.390e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.8870e+00 J
  → Fracture energy : 1.5712e+00 J
  → Total energy    : 6.4583e+00 J


## Step 77/151: t = 1.82e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 76 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.52e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.451e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.089e-02
  |ΔD|_∞ = 1.101e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.52e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.177e-03
  |ΔD|_∞ = 2.202e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.52e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.355e-04
  |ΔD|_∞ = 4.404e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.52e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.710e-05
  |ΔD|_∞ = 8.807e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.52e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.742e-05
  |ΔD|_∞ = 1.761e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.52e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.484e-06
  |ΔD|_∞ = 3.523e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.9933e+00 J
  → Fracture energy : 1.5844e+00 J
  → Total energy    : 6.5777e+00 J


## Step 78/151: t = 1.85e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 77 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.54e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.448e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.089e-02
  |ΔD|_∞ = 1.133e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.54e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.177e-03
  |ΔD|_∞ = 2.266e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.54e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.354e-04
  |ΔD|_∞ = 4.532e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.54e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.708e-05
  |ΔD|_∞ = 9.065e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.54e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.742e-05
  |ΔD|_∞ = 1.813e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.54e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.483e-06
  |ΔD|_∞ = 3.626e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.0997e+00 J
  → Fracture energy : 1.5984e+00 J
  → Total energy    : 6.6981e+00 J


## Step 79/151: t = 1.87e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 78 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.56e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.429e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.088e-02
  |ΔD|_∞ = 1.142e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.56e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.176e-03
  |ΔD|_∞ = 2.284e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.56e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.353e-04
  |ΔD|_∞ = 4.567e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.56e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.705e-05
  |ΔD|_∞ = 9.135e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.56e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.741e-05
  |ΔD|_∞ = 1.827e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.56e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.482e-06
  |ΔD|_∞ = 3.654e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.2044e+00 J
  → Fracture energy : 1.6131e+00 J
  → Total energy    : 6.8175e+00 J


## Step 80/151: t = 1.90e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 79 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.58e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.463e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.083e-02
  |ΔD|_∞ = 1.173e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.58e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.166e-03
  |ΔD|_∞ = 2.345e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.58e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.206e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.332e-04
  |ΔD|_∞ = 4.691e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.58e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.206e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.664e-05
  |ΔD|_∞ = 9.382e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.58e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.733e-05
  |ΔD|_∞ = 1.876e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.58e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.465e-06
  |ΔD|_∞ = 3.753e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.3086e+00 J
  → Fracture energy : 1.6285e+00 J
  → Total energy    : 6.9372e+00 J


## Step 81/151: t = 1.92e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 80 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.6e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.428e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.069e-02
  |ΔD|_∞ = 1.182e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.6e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.139e-03
  |ΔD|_∞ = 2.364e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.6e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.277e-04
  |ΔD|_∞ = 4.728e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.6e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.555e-05
  |ΔD|_∞ = 9.457e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.6e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.711e-05
  |ΔD|_∞ = 1.891e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.6e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.386e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.422e-06
  |ΔD|_∞ = 3.783e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.4126e+00 J
  → Fracture energy : 1.6447e+00 J
  → Total energy    : 7.0573e+00 J


## Step 82/151: t = 1.94e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 81 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.62e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.428e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.047e-02
  |ΔD|_∞ = 1.189e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.62e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.093e-03
  |ΔD|_∞ = 2.378e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.62e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.187e-04
  |ΔD|_∞ = 4.755e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.62e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.373e-05
  |ΔD|_∞ = 9.510e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.62e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.675e-05
  |ΔD|_∞ = 1.902e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.62e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.349e-06
  |ΔD|_∞ = 3.804e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.5169e+00 J
  → Fracture energy : 1.6616e+00 J
  → Total energy    : 7.1784e+00 J


## Step 83/151: t = 1.97e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 82 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.64e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.380e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.033e-02
  |ΔD|_∞ = 1.202e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.64e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.836e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.065e-03
  |ΔD|_∞ = 2.405e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.64e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.836e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.130e-04
  |ΔD|_∞ = 4.810e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.64e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.260e-05
  |ΔD|_∞ = 9.620e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.64e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.652e-05
  |ΔD|_∞ = 1.924e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.64e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.304e-06
  |ΔD|_∞ = 3.848e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.6213e+00 J
  → Fracture energy : 1.6791e+00 J
  → Total energy    : 7.3004e+00 J


## Step 84/151: t = 1.99e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 83 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.66e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.385e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.027e-02
  |ΔD|_∞ = 1.215e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.66e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.054e-03
  |ΔD|_∞ = 2.431e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.66e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.108e-04
  |ΔD|_∞ = 4.862e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.66e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.216e-05
  |ΔD|_∞ = 9.723e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.66e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.643e-05
  |ΔD|_∞ = 1.945e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.66e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.407e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.286e-06
  |ΔD|_∞ = 3.889e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.7240e+00 J
  → Fracture energy : 1.6975e+00 J
  → Total energy    : 7.4215e+00 J


## Step 85/151: t = 2.02e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 84 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.68e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.376e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.020e-02
  |ΔD|_∞ = 1.250e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.68e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.233e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.040e-03
  |ΔD|_∞ = 2.499e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.68e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.233e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.080e-04
  |ΔD|_∞ = 4.998e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.68e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.160e-05
  |ΔD|_∞ = 9.996e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.68e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.632e-05
  |ΔD|_∞ = 1.999e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.68e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.264e-06
  |ΔD|_∞ = 3.999e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.8273e+00 J
  → Fracture energy : 1.7165e+00 J
  → Total energy    : 7.5439e+00 J


## Step 86/151: t = 2.04e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 85 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.7e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.353e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.007e-02
  |ΔD|_∞ = 1.264e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.7e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.013e-03
  |ΔD|_∞ = 2.528e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.7e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.027e-04
  |ΔD|_∞ = 5.056e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.7e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.053e-05
  |ΔD|_∞ = 1.011e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.7e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.611e-05
  |ΔD|_∞ = 2.022e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.7e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.671e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.221e-06
  |ΔD|_∞ = 4.045e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.9290e+00 J
  → Fracture energy : 1.7364e+00 J
  → Total energy    : 7.6654e+00 J


## Step 87/151: t = 2.06e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 86 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.72e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.355e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.917e-03
  |ΔD|_∞ = 1.275e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.72e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.983e-03
  |ΔD|_∞ = 2.550e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.72e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.967e-04
  |ΔD|_∞ = 5.099e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.72e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.934e-05
  |ΔD|_∞ = 1.020e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.72e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.587e-05
  |ΔD|_∞ = 2.040e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.72e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.174e-06
  |ΔD|_∞ = 4.079e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.0281e+00 J
  → Fracture energy : 1.7571e+00 J
  → Total energy    : 7.7851e+00 J


## Step 88/151: t = 2.09e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 87 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.74e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.378e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.841e-03
  |ΔD|_∞ = 1.301e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.74e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.198e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.968e-03
  |ΔD|_∞ = 2.601e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.74e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.936e-04
  |ΔD|_∞ = 5.202e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.74e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.872e-05
  |ΔD|_∞ = 1.040e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.74e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.574e-05
  |ΔD|_∞ = 2.081e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.74e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.149e-06
  |ΔD|_∞ = 4.162e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.1280e+00 J
  → Fracture energy : 1.7786e+00 J
  → Total energy    : 7.9066e+00 J


## Step 89/151: t = 2.11e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 88 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.76e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.331e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.798e-03
  |ΔD|_∞ = 1.307e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.76e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.960e-03
  |ΔD|_∞ = 2.615e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.76e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.919e-04
  |ΔD|_∞ = 5.230e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.76e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.838e-05
  |ΔD|_∞ = 1.046e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.76e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.669e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.568e-05
  |ΔD|_∞ = 2.092e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.76e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.519e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.135e-06
  |ΔD|_∞ = 4.184e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.2288e+00 J
  → Fracture energy : 1.8009e+00 J
  → Total energy    : 8.0297e+00 J


## Step 90/151: t = 2.14e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 89 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.78e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.313e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.765e-03
  |ΔD|_∞ = 1.298e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.78e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.953e-03
  |ΔD|_∞ = 2.595e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.78e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.650e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.906e-04
  |ΔD|_∞ = 5.190e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.78e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.650e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.812e-05
  |ΔD|_∞ = 1.038e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.78e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.562e-05
  |ΔD|_∞ = 2.076e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.78e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.125e-06
  |ΔD|_∞ = 4.152e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.3288e+00 J
  → Fracture energy : 1.8240e+00 J
  → Total energy    : 8.1529e+00 J


## Step 91/151: t = 2.16e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 90 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.8e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.307e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.688e-03
  |ΔD|_∞ = 1.324e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.8e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.938e-03
  |ΔD|_∞ = 2.648e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.8e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.875e-04
  |ΔD|_∞ = 5.297e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.8e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.750e-05
  |ΔD|_∞ = 1.059e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.8e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.550e-05
  |ΔD|_∞ = 2.119e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.8e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.100e-06
  |ΔD|_∞ = 4.237e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.4281e+00 J
  → Fracture energy : 1.8481e+00 J
  → Total energy    : 8.2762e+00 J


## Step 92/151: t = 2.18e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 91 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.82e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.297e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.597e-03
  |ΔD|_∞ = 1.353e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.82e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.919e-03
  |ΔD|_∞ = 2.706e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.82e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.839e-04
  |ΔD|_∞ = 5.412e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.82e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.678e-05
  |ΔD|_∞ = 1.082e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.82e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.536e-05
  |ΔD|_∞ = 2.165e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.82e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.501e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.071e-06
  |ΔD|_∞ = 4.330e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.5265e+00 J
  → Fracture energy : 1.8729e+00 J
  → Total energy    : 8.3994e+00 J


## Step 93/151: t = 2.21e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 92 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.84e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.290e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.598e-03
  |ΔD|_∞ = 1.338e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.84e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.920e-03
  |ΔD|_∞ = 2.677e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.84e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.839e-04
  |ΔD|_∞ = 5.354e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.84e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.678e-05
  |ΔD|_∞ = 1.071e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.84e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.536e-05
  |ΔD|_∞ = 2.142e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.84e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.071e-06
  |ΔD|_∞ = 4.283e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.6235e+00 J
  → Fracture energy : 1.8984e+00 J
  → Total energy    : 8.5219e+00 J


## Step 94/151: t = 2.23e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 93 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.86e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.289e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.670e-03
  |ΔD|_∞ = 1.352e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.86e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.039e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.934e-03
  |ΔD|_∞ = 2.705e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.86e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.056e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.868e-04
  |ΔD|_∞ = 5.410e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.86e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.736e-05
  |ΔD|_∞ = 1.082e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.86e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.547e-05
  |ΔD|_∞ = 2.164e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.86e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.095e-06
  |ΔD|_∞ = 4.328e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.7190e+00 J
  → Fracture energy : 1.9248e+00 J
  → Total energy    : 8.6438e+00 J


## Step 95/151: t = 2.26e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 94 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.88e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.282e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.678e-03
  |ΔD|_∞ = 1.382e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.88e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.643e-20

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.936e-03
  |ΔD|_∞ = 2.764e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.88e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.871e-04
  |ΔD|_∞ = 5.527e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.88e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.742e-05
  |ΔD|_∞ = 1.105e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.88e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.548e-05
  |ΔD|_∞ = 2.211e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.88e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.188e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.097e-06
  |ΔD|_∞ = 4.422e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.8131e+00 J
  → Fracture energy : 1.9523e+00 J
  → Total energy    : 8.7654e+00 J


## Step 96/151: t = 2.28e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 95 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.9e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.280e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.544e-03
  |ΔD|_∞ = 1.361e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.9e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.909e-03
  |ΔD|_∞ = 2.721e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.9e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.818e-04
  |ΔD|_∞ = 5.443e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.9e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.635e-05
  |ΔD|_∞ = 1.089e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.9e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.527e-05
  |ΔD|_∞ = 2.177e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.9e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.054e-06
  |ΔD|_∞ = 4.354e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.9066e+00 J
  → Fracture energy : 1.9805e+00 J
  → Total energy    : 8.8871e+00 J


## Step 97/151: t = 2.30e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 96 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.92e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.268e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.434e-03
  |ΔD|_∞ = 1.381e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.92e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.887e-03
  |ΔD|_∞ = 2.762e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.92e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.774e-04
  |ΔD|_∞ = 5.525e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.92e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.569e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.547e-05
  |ΔD|_∞ = 1.105e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.92e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.569e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.509e-05
  |ΔD|_∞ = 2.210e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.92e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.019e-06
  |ΔD|_∞ = 4.420e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.9987e+00 J
  → Fracture energy : 2.0096e+00 J
  → Total energy    : 9.0084e+00 J


## Step 98/151: t = 2.33e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 97 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.94e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.262e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.383e-03
  |ΔD|_∞ = 1.402e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.94e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.877e-03
  |ΔD|_∞ = 2.804e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.94e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.753e-04
  |ΔD|_∞ = 5.609e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.94e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.507e-05
  |ΔD|_∞ = 1.122e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.94e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.501e-05
  |ΔD|_∞ = 2.243e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.94e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.448e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.003e-06
  |ΔD|_∞ = 4.487e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.0895e+00 J
  → Fracture energy : 2.0395e+00 J
  → Total energy    : 9.1291e+00 J


## Step 99/151: t = 2.35e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 98 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.96e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.257e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.317e-03
  |ΔD|_∞ = 1.416e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.96e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.863e-03
  |ΔD|_∞ = 2.832e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.96e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.727e-04
  |ΔD|_∞ = 5.663e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.96e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.454e-05
  |ΔD|_∞ = 1.133e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.96e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.491e-05
  |ΔD|_∞ = 2.265e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.96e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.982e-06
  |ΔD|_∞ = 4.530e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.1788e+00 J
  → Fracture energy : 2.0704e+00 J
  → Total energy    : 9.2492e+00 J


## Step 100/151: t = 2.38e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 99 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.98e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.251e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.131e-03
  |ΔD|_∞ = 1.390e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.98e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.826e-03
  |ΔD|_∞ = 2.779e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.98e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.652e-04
  |ΔD|_∞ = 5.558e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.98e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.305e-05
  |ΔD|_∞ = 1.112e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.98e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.461e-05
  |ΔD|_∞ = 2.223e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [1.98e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.922e-06
  |ΔD|_∞ = 4.447e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.2668e+00 J
  → Fracture energy : 2.1021e+00 J
  → Total energy    : 9.3690e+00 J


## Step 101/151: t = 2.40e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 100 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.245e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.933e-03
  |ΔD|_∞ = 1.399e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.787e-03
  |ΔD|_∞ = 2.799e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.573e-04
  |ΔD|_∞ = 5.598e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.146e-05
  |ΔD|_∞ = 1.120e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.429e-05
  |ΔD|_∞ = 2.239e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.858e-06
  |ΔD|_∞ = 4.478e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.3534e+00 J
  → Fracture energy : 2.1348e+00 J
  → Total energy    : 9.4882e+00 J


## Step 102/151: t = 2.42e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 101 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.02e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.238e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.855e-03
  |ΔD|_∞ = 1.399e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.02e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.045e-20

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.771e-03
  |ΔD|_∞ = 2.798e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.02e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.542e-04
  |ΔD|_∞ = 5.596e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.02e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.084e-05
  |ΔD|_∞ = 1.119e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.02e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.417e-05
  |ΔD|_∞ = 2.238e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.02e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.833e-06
  |ΔD|_∞ = 4.477e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.4386e+00 J
  → Fracture energy : 2.1684e+00 J
  → Total energy    : 9.6070e+00 J


## Step 103/151: t = 2.45e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 102 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.04e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.234e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.865e-03
  |ΔD|_∞ = 1.396e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.04e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.773e-03
  |ΔD|_∞ = 2.793e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.04e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.546e-04
  |ΔD|_∞ = 5.586e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.04e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.092e-05
  |ΔD|_∞ = 1.117e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.04e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.418e-05
  |ΔD|_∞ = 2.234e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.04e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.837e-06
  |ΔD|_∞ = 4.469e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.5222e+00 J
  → Fracture energy : 2.2027e+00 J
  → Total energy    : 9.7249e+00 J


## Step 104/151: t = 2.47e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 103 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.06e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.231e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.934e-03
  |ΔD|_∞ = 1.401e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.06e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.787e-03
  |ΔD|_∞ = 2.802e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.06e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.574e-04
  |ΔD|_∞ = 5.605e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.06e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.148e-05
  |ΔD|_∞ = 1.121e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.06e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.430e-05
  |ΔD|_∞ = 2.242e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.06e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.859e-06
  |ΔD|_∞ = 4.484e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.6036e+00 J
  → Fracture energy : 2.2379e+00 J
  → Total energy    : 9.8416e+00 J


## Step 105/151: t = 2.50e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 104 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.08e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.225e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.961e-03
  |ΔD|_∞ = 1.412e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.08e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.792e-03
  |ΔD|_∞ = 2.823e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.08e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.584e-04
  |ΔD|_∞ = 5.646e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.08e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.169e-05
  |ΔD|_∞ = 1.129e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.08e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.434e-05
  |ΔD|_∞ = 2.259e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.08e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.867e-06
  |ΔD|_∞ = 4.517e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.6835e+00 J
  → Fracture energy : 2.2740e+00 J
  → Total energy    : 9.9575e+00 J


## Step 106/151: t = 2.52e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 105 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.1e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.219e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.894e-03
  |ΔD|_∞ = 1.405e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.1e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.422e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.779e-03
  |ΔD|_∞ = 2.811e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.1e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.557e-04
  |ΔD|_∞ = 5.621e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.1e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.468e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.115e-05
  |ΔD|_∞ = 1.124e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.1e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.468e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.423e-05
  |ΔD|_∞ = 2.249e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.1e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.846e-06
  |ΔD|_∞ = 4.497e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.7615e+00 J
  → Fracture energy : 2.3110e+00 J
  → Total energy    : 1.0072e+01 J


## Step 107/151: t = 2.54e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 106 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.12e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.214e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.713e-03
  |ΔD|_∞ = 1.418e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.12e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.743e-03
  |ΔD|_∞ = 2.836e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.12e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.485e-04
  |ΔD|_∞ = 5.673e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.12e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.971e-05
  |ΔD|_∞ = 1.135e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.12e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.394e-05
  |ΔD|_∞ = 2.269e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.12e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.788e-06
  |ΔD|_∞ = 4.538e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.8378e+00 J
  → Fracture energy : 2.3489e+00 J
  → Total energy    : 1.0187e+01 J


## Step 108/151: t = 2.57e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 107 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.14e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.209e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.534e-03
  |ΔD|_∞ = 1.466e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.14e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.707e-03
  |ΔD|_∞ = 2.931e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.14e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.414e-04
  |ΔD|_∞ = 5.862e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.14e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.828e-05
  |ΔD|_∞ = 1.172e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.14e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.555e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.366e-05
  |ΔD|_∞ = 2.345e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.14e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.555e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.731e-06
  |ΔD|_∞ = 4.690e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.9122e+00 J
  → Fracture energy : 2.3877e+00 J
  → Total energy    : 1.0300e+01 J


## Step 109/151: t = 2.59e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 108 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.16e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.206e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.498e-03
  |ΔD|_∞ = 1.448e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.16e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.700e-03
  |ΔD|_∞ = 2.897e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.16e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.218e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.399e-04
  |ΔD|_∞ = 5.794e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.16e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.589e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.799e-05
  |ΔD|_∞ = 1.159e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.16e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.500e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.360e-05
  |ΔD|_∞ = 2.317e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.16e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.221e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.719e-06
  |ΔD|_∞ = 4.635e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.9843e+00 J
  → Fracture energy : 2.4273e+00 J
  → Total energy    : 1.0412e+01 J


## Step 110/151: t = 2.62e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 109 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.18e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.198e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.575e-03
  |ΔD|_∞ = 1.448e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.18e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.715e-03
  |ΔD|_∞ = 2.897e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.18e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.430e-04
  |ΔD|_∞ = 5.794e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.18e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.860e-05
  |ΔD|_∞ = 1.159e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.18e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.372e-05
  |ΔD|_∞ = 2.317e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.18e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.744e-06
  |ΔD|_∞ = 4.635e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.0548e+00 J
  → Fracture energy : 2.4683e+00 J
  → Total energy    : 1.0523e+01 J


## Step 111/151: t = 2.64e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 110 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.2e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.190e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.601e-03
  |ΔD|_∞ = 1.439e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.2e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.720e-03
  |ΔD|_∞ = 2.878e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.2e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.440e-04
  |ΔD|_∞ = 5.756e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.2e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.096e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.881e-05
  |ΔD|_∞ = 1.151e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.2e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.096e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.376e-05
  |ΔD|_∞ = 2.302e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.2e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.752e-06
  |ΔD|_∞ = 4.604e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.1237e+00 J
  → Fracture energy : 2.5104e+00 J
  → Total energy    : 1.0634e+01 J


## Step 112/151: t = 2.66e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 111 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.22e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.183e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.492e-03
  |ΔD|_∞ = 1.423e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.22e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.698e-03
  |ΔD|_∞ = 2.845e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.22e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.397e-04
  |ΔD|_∞ = 5.690e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.22e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.794e-05
  |ΔD|_∞ = 1.138e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.22e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.359e-05
  |ΔD|_∞ = 2.276e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.22e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.717e-06
  |ΔD|_∞ = 4.552e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.1908e+00 J
  → Fracture energy : 2.5533e+00 J
  → Total energy    : 1.0744e+01 J


## Step 113/151: t = 2.69e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 112 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.24e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.181e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.249e-03
  |ΔD|_∞ = 1.384e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.24e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.650e-03
  |ΔD|_∞ = 2.768e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.24e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.300e-04
  |ΔD|_∞ = 5.536e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.24e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.138e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.599e-05
  |ΔD|_∞ = 1.107e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.24e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.778e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.320e-05
  |ΔD|_∞ = 2.214e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.24e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.791e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.640e-06
  |ΔD|_∞ = 4.429e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.2567e+00 J
  → Fracture energy : 2.5969e+00 J
  → Total energy    : 1.0854e+01 J


## Step 114/151: t = 2.71e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 113 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.26e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.166e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.037e-03
  |ΔD|_∞ = 1.415e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.26e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.607e-03
  |ΔD|_∞ = 2.831e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.26e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.215e-04
  |ΔD|_∞ = 5.662e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.26e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.429e-05
  |ΔD|_∞ = 1.132e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.26e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.286e-05
  |ΔD|_∞ = 2.265e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.26e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.488e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.572e-06
  |ΔD|_∞ = 4.529e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.3206e+00 J
  → Fracture energy : 2.6415e+00 J
  → Total energy    : 1.0962e+01 J


## Step 115/151: t = 2.74e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 114 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.28e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.161e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.989e-03
  |ΔD|_∞ = 1.428e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.28e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.598e-03
  |ΔD|_∞ = 2.855e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.28e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.196e-04
  |ΔD|_∞ = 5.711e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.28e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.391e-05
  |ΔD|_∞ = 1.142e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.28e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.278e-05
  |ΔD|_∞ = 2.284e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.28e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.556e-06
  |ΔD|_∞ = 4.569e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.3822e+00 J
  → Fracture energy : 2.6871e+00 J
  → Total energy    : 1.1069e+01 J


## Step 116/151: t = 2.76e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 115 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.3e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.160e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.062e-03
  |ΔD|_∞ = 1.446e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.3e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.645e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.612e-03
  |ΔD|_∞ = 2.892e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.3e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.225e-04
  |ΔD|_∞ = 5.785e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.3e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.450e-05
  |ΔD|_∞ = 1.157e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.3e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.290e-05
  |ΔD|_∞ = 2.314e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.3e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.580e-06
  |ΔD|_∞ = 4.628e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.4385e+00 J
  → Fracture energy : 2.7333e+00 J
  → Total energy    : 1.1172e+01 J


## Step 117/151: t = 2.78e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 116 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.32e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.189e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.176e-03
  |ΔD|_∞ = 1.460e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.32e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.635e-03
  |ΔD|_∞ = 2.920e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.32e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.270e-04
  |ΔD|_∞ = 5.840e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.32e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.541e-05
  |ΔD|_∞ = 1.168e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.32e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.308e-05
  |ΔD|_∞ = 2.336e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.32e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.616e-06
  |ΔD|_∞ = 4.672e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.4912e+00 J
  → Fracture energy : 2.7804e+00 J
  → Total energy    : 1.1272e+01 J


## Step 118/151: t = 2.81e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 117 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.34e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.178e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.234e-03
  |ΔD|_∞ = 1.463e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.34e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.647e-03
  |ΔD|_∞ = 2.926e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.34e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.294e-04
  |ΔD|_∞ = 5.852e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.34e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.049e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.587e-05
  |ΔD|_∞ = 1.170e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.34e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.049e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.317e-05
  |ΔD|_∞ = 2.341e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.34e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.635e-06
  |ΔD|_∞ = 4.681e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.5453e+00 J
  → Fracture energy : 2.8281e+00 J
  → Total energy    : 1.1373e+01 J


## Step 119/151: t = 2.83e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 118 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.36e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.139e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.282e-03
  |ΔD|_∞ = 1.492e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.36e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.297e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.656e-03
  |ΔD|_∞ = 2.984e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.36e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.313e-04
  |ΔD|_∞ = 5.969e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.36e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.625e-05
  |ΔD|_∞ = 1.194e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.36e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.325e-05
  |ΔD|_∞ = 2.387e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.36e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.650e-06
  |ΔD|_∞ = 4.775e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.5979e+00 J
  → Fracture energy : 2.8768e+00 J
  → Total energy    : 1.1475e+01 J


## Step 120/151: t = 2.86e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 119 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.38e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.131e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.271e-03
  |ΔD|_∞ = 1.489e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.38e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.654e-03
  |ΔD|_∞ = 2.978e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.38e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.309e-04
  |ΔD|_∞ = 5.956e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.38e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.617e-05
  |ΔD|_∞ = 1.191e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.38e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.323e-05
  |ΔD|_∞ = 2.382e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.38e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.647e-06
  |ΔD|_∞ = 4.765e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.6488e+00 J
  → Fracture energy : 2.9262e+00 J
  → Total energy    : 1.1575e+01 J


## Step 121/151: t = 2.88e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 120 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.4e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.125e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.182e-03
  |ΔD|_∞ = 1.469e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.4e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.636e-03
  |ΔD|_∞ = 2.938e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.4e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.273e-04
  |ΔD|_∞ = 5.877e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.4e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.545e-05
  |ΔD|_∞ = 1.175e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.4e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.951e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.309e-05
  |ΔD|_∞ = 2.351e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.4e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.951e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.618e-06
  |ΔD|_∞ = 4.701e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.6966e+00 J
  → Fracture energy : 2.9770e+00 J
  → Total energy    : 1.1674e+01 J


## Step 122/151: t = 2.90e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 121 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.42e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.119e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.975e-03
  |ΔD|_∞ = 1.476e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.42e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.595e-03
  |ΔD|_∞ = 2.952e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.42e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.861e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.190e-04
  |ΔD|_∞ = 5.905e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.42e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.861e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.380e-05
  |ΔD|_∞ = 1.181e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.42e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.276e-05
  |ΔD|_∞ = 2.362e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.42e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.552e-06
  |ΔD|_∞ = 4.724e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.7428e+00 J
  → Fracture energy : 3.0287e+00 J
  → Total energy    : 1.1772e+01 J


## Step 123/151: t = 2.93e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 122 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.44e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.109e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.658e-03
  |ΔD|_∞ = 1.453e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.44e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.532e-03
  |ΔD|_∞ = 2.906e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.44e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.063e-04
  |ΔD|_∞ = 5.813e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.44e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.127e-05
  |ΔD|_∞ = 1.163e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.44e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.225e-05
  |ΔD|_∞ = 2.325e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.44e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.451e-06
  |ΔD|_∞ = 4.650e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.7886e+00 J
  → Fracture energy : 3.0804e+00 J
  → Total energy    : 1.1869e+01 J


## Step 124/151: t = 2.95e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 123 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.46e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.101e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.530e-03
  |ΔD|_∞ = 1.435e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.46e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.079e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.506e-03
  |ΔD|_∞ = 2.869e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.46e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.079e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.012e-04
  |ΔD|_∞ = 5.739e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.46e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.024e-05
  |ΔD|_∞ = 1.148e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.46e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.205e-05
  |ΔD|_∞ = 2.296e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.46e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.410e-06
  |ΔD|_∞ = 4.591e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.8323e+00 J
  → Fracture energy : 3.1327e+00 J
  → Total energy    : 1.1965e+01 J


## Step 125/151: t = 2.98e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 124 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.48e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.093e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.527e-03
  |ΔD|_∞ = 1.476e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.48e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.505e-03
  |ΔD|_∞ = 2.952e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.48e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.011e-04
  |ΔD|_∞ = 5.905e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.48e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.021e-05
  |ΔD|_∞ = 1.181e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.48e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.204e-05
  |ΔD|_∞ = 2.362e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.48e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.409e-06
  |ΔD|_∞ = 4.724e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.8743e+00 J
  → Fracture energy : 3.1856e+00 J
  → Total energy    : 1.2060e+01 J


## Step 126/151: t = 3.00e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 125 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.5e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.085e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.560e-03
  |ΔD|_∞ = 1.475e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.5e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.512e-03
  |ΔD|_∞ = 2.950e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.5e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.024e-04
  |ΔD|_∞ = 5.901e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.5e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.048e-05
  |ΔD|_∞ = 1.180e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.5e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.210e-05
  |ΔD|_∞ = 2.360e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.5e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.419e-06
  |ΔD|_∞ = 4.720e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.9124e+00 J
  → Fracture energy : 3.2391e+00 J
  → Total energy    : 1.2151e+01 J


## Step 127/151: t = 3.02e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 126 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.52e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.101e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.542e-03
  |ΔD|_∞ = 1.475e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.52e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.945e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.508e-03
  |ΔD|_∞ = 2.950e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.52e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.017e-04
  |ΔD|_∞ = 5.900e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.52e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.034e-05
  |ΔD|_∞ = 1.180e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.52e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.207e-05
  |ΔD|_∞ = 2.360e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.52e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.414e-06
  |ΔD|_∞ = 4.720e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.9501e+00 J
  → Fracture energy : 3.2932e+00 J
  → Total energy    : 1.2243e+01 J


## Step 128/151: t = 3.05e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 127 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.54e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.067e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.499e-03
  |ΔD|_∞ = 1.468e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.54e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.619e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.500e-03
  |ΔD|_∞ = 2.936e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.54e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.000e-04
  |ΔD|_∞ = 5.871e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.54e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.999e-05
  |ΔD|_∞ = 1.174e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.54e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.200e-05
  |ΔD|_∞ = 2.349e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.54e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.400e-06
  |ΔD|_∞ = 4.697e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.9873e+00 J
  → Fracture energy : 3.3480e+00 J
  → Total energy    : 1.2335e+01 J


## Step 129/151: t = 3.07e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 128 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.56e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.058e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.403e-03
  |ΔD|_∞ = 1.446e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.56e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.853e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.481e-03
  |ΔD|_∞ = 2.892e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.56e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.853e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.961e-04
  |ΔD|_∞ = 5.784e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.56e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.912e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.922e-05
  |ΔD|_∞ = 1.157e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.56e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.912e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.184e-05
  |ΔD|_∞ = 2.313e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.56e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.369e-06
  |ΔD|_∞ = 4.627e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.0238e+00 J
  → Fracture energy : 3.4033e+00 J
  → Total energy    : 1.2427e+01 J


## Step 130/151: t = 3.10e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 129 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.58e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.048e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.278e-03
  |ΔD|_∞ = 1.442e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.58e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.456e-03
  |ΔD|_∞ = 2.884e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.58e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.911e-04
  |ΔD|_∞ = 5.768e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.58e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.823e-05
  |ΔD|_∞ = 1.154e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.58e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.165e-05
  |ΔD|_∞ = 2.307e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.58e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.329e-06
  |ΔD|_∞ = 4.614e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.0589e+00 J
  → Fracture energy : 3.4597e+00 J
  → Total energy    : 1.2519e+01 J


## Step 131/151: t = 3.12e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 130 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.6e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.039e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.166e-03
  |ΔD|_∞ = 1.462e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.6e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.352e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.433e-03
  |ΔD|_∞ = 2.925e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.6e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.352e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.866e-04
  |ΔD|_∞ = 5.850e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.6e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.733e-05
  |ΔD|_∞ = 1.170e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.6e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.147e-05
  |ΔD|_∞ = 2.340e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.6e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.293e-06
  |ΔD|_∞ = 4.680e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.0930e+00 J
  → Fracture energy : 3.5165e+00 J
  → Total energy    : 1.2609e+01 J


## Step 132/151: t = 3.14e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 131 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.62e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.033e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.094e-03
  |ΔD|_∞ = 1.461e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.62e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.419e-03
  |ΔD|_∞ = 2.923e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.62e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.837e-04
  |ΔD|_∞ = 5.845e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.62e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.675e-05
  |ΔD|_∞ = 1.169e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.62e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.135e-05
  |ΔD|_∞ = 2.338e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.62e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.270e-06
  |ΔD|_∞ = 4.676e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.1255e+00 J
  → Fracture energy : 3.5734e+00 J
  → Total energy    : 1.2699e+01 J


## Step 133/151: t = 3.17e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 132 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.64e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.025e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.061e-03
  |ΔD|_∞ = 1.454e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.64e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.412e-03
  |ΔD|_∞ = 2.908e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.64e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.824e-04
  |ΔD|_∞ = 5.815e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.64e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.649e-05
  |ΔD|_∞ = 1.163e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.64e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.130e-05
  |ΔD|_∞ = 2.326e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.64e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.259e-06
  |ΔD|_∞ = 4.652e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.1564e+00 J
  → Fracture energy : 3.6306e+00 J
  → Total energy    : 1.2787e+01 J


## Step 134/151: t = 3.19e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 133 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.66e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.015e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.010e-03
  |ΔD|_∞ = 1.431e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.66e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.402e-03
  |ΔD|_∞ = 2.863e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.66e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.804e-04
  |ΔD|_∞ = 5.726e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.66e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.608e-05
  |ΔD|_∞ = 1.145e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.66e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.122e-05
  |ΔD|_∞ = 2.290e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.66e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.243e-06
  |ΔD|_∞ = 4.580e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.1864e+00 J
  → Fracture energy : 3.6881e+00 J
  → Total energy    : 1.2874e+01 J


## Step 135/151: t = 3.22e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 134 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.68e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.005e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.938e-03
  |ΔD|_∞ = 1.479e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.68e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.930e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.388e-03
  |ΔD|_∞ = 2.959e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.68e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.775e-04
  |ΔD|_∞ = 5.917e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.68e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.551e-05
  |ΔD|_∞ = 1.183e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.68e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.110e-05
  |ΔD|_∞ = 2.367e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.68e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.220e-06
  |ΔD|_∞ = 4.734e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.2157e+00 J
  → Fracture energy : 3.7464e+00 J
  → Total energy    : 1.2962e+01 J


## Step 136/151: t = 3.24e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 135 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.7e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.973e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.851e-03
  |ΔD|_∞ = 1.496e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.7e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.226e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.370e-03
  |ΔD|_∞ = 2.993e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.7e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.740e-04
  |ΔD|_∞ = 5.986e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.7e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.481e-05
  |ΔD|_∞ = 1.197e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.7e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.224e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.096e-05
  |ΔD|_∞ = 2.394e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.7e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.224e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.192e-06
  |ΔD|_∞ = 4.789e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.2436e+00 J
  → Fracture energy : 3.8049e+00 J
  → Total energy    : 1.3048e+01 J


## Step 137/151: t = 3.26e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 136 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.72e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.918e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.761e-03
  |ΔD|_∞ = 1.485e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.72e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.352e-03
  |ΔD|_∞ = 2.971e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.72e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.704e-04
  |ΔD|_∞ = 5.942e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.72e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.409e-05
  |ΔD|_∞ = 1.188e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.72e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.082e-05
  |ΔD|_∞ = 2.377e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.72e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.164e-06
  |ΔD|_∞ = 4.753e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.2702e+00 J
  → Fracture energy : 3.8628e+00 J
  → Total energy    : 1.3133e+01 J


## Step 138/151: t = 3.29e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 137 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.74e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.827e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.693e-03
  |ΔD|_∞ = 1.430e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.74e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.339e-03
  |ΔD|_∞ = 2.860e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.74e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.677e-04
  |ΔD|_∞ = 5.719e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.74e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.354e-05
  |ΔD|_∞ = 1.144e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.74e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.071e-05
  |ΔD|_∞ = 2.288e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.74e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.142e-06
  |ΔD|_∞ = 4.575e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.2957e+00 J
  → Fracture energy : 3.9210e+00 J
  → Total energy    : 1.3217e+01 J


## Step 139/151: t = 3.31e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 138 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.76e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.722e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.636e-03
  |ΔD|_∞ = 1.474e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.76e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.327e-03
  |ΔD|_∞ = 2.948e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.76e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.654e-04
  |ΔD|_∞ = 5.896e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.76e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.309e-05
  |ΔD|_∞ = 1.179e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.76e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.062e-05
  |ΔD|_∞ = 2.358e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.76e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.753e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.123e-06
  |ΔD|_∞ = 4.717e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.3206e+00 J
  → Fracture energy : 3.9798e+00 J
  → Total energy    : 1.3300e+01 J


## Step 140/151: t = 3.34e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 139 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.78e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.637e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.651e-03
  |ΔD|_∞ = 1.509e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.78e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.330e-03
  |ΔD|_∞ = 3.019e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.78e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.661e-04
  |ΔD|_∞ = 6.037e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.78e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.321e-05
  |ΔD|_∞ = 1.207e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.78e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.064e-05
  |ΔD|_∞ = 2.415e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.78e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.128e-06
  |ΔD|_∞ = 4.830e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.3444e+00 J
  → Fracture energy : 4.0400e+00 J
  → Total energy    : 1.3384e+01 J


## Step 141/151: t = 3.36e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 140 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.8e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.574e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.577e-03
  |ΔD|_∞ = 1.435e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.8e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.315e-03
  |ΔD|_∞ = 2.870e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.8e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.631e-04
  |ΔD|_∞ = 5.741e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.8e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.261e-05
  |ΔD|_∞ = 1.148e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.8e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.052e-05
  |ΔD|_∞ = 2.296e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.8e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.105e-06
  |ΔD|_∞ = 4.593e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.3697e+00 J
  → Fracture energy : 4.0981e+00 J
  → Total energy    : 1.3468e+01 J


## Step 142/151: t = 3.38e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 141 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.82e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.478e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.529e-03
  |ΔD|_∞ = 1.419e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.82e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.306e-03
  |ΔD|_∞ = 2.837e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.82e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.612e-04
  |ΔD|_∞ = 5.675e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.82e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.223e-05
  |ΔD|_∞ = 1.135e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.82e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.045e-05
  |ΔD|_∞ = 2.270e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.82e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.853e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.089e-06
  |ΔD|_∞ = 4.540e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.3906e+00 J
  → Fracture energy : 4.1584e+00 J
  → Total energy    : 1.3549e+01 J


## Step 143/151: t = 3.41e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 142 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.84e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.433e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.437e-03
  |ΔD|_∞ = 1.469e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.84e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.287e-03
  |ΔD|_∞ = 2.938e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.84e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.575e-04
  |ΔD|_∞ = 5.877e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.84e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.149e-05
  |ΔD|_∞ = 1.175e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.84e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.512e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.030e-05
  |ΔD|_∞ = 2.351e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.84e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.512e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.060e-06
  |ΔD|_∞ = 4.702e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.4128e+00 J
  → Fracture energy : 4.2179e+00 J
  → Total energy    : 1.3631e+01 J


## Step 144/151: t = 3.43e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 143 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.86e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.369e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.470e-03
  |ΔD|_∞ = 1.451e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.86e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.971e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.294e-03
  |ΔD|_∞ = 2.902e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.86e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.971e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.588e-04
  |ΔD|_∞ = 5.805e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.86e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.176e-05
  |ΔD|_∞ = 1.161e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.86e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.035e-05
  |ΔD|_∞ = 2.322e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.86e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.070e-06
  |ΔD|_∞ = 4.644e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.4302e+00 J
  → Fracture energy : 4.2787e+00 J
  → Total energy    : 1.3709e+01 J


## Step 145/151: t = 3.46e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 144 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.88e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.339e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.412e-03
  |ΔD|_∞ = 1.489e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.88e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.282e-03
  |ΔD|_∞ = 2.979e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.88e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.565e-04
  |ΔD|_∞ = 5.958e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.88e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.130e-05
  |ΔD|_∞ = 1.192e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.88e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.026e-05
  |ΔD|_∞ = 2.383e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.88e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.052e-06
  |ΔD|_∞ = 4.766e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.4495e+00 J
  → Fracture energy : 4.3375e+00 J
  → Total energy    : 1.3787e+01 J


## Step 146/151: t = 3.48e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 145 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.9e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.238e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.383e-03
  |ΔD|_∞ = 1.522e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.9e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.277e-03
  |ΔD|_∞ = 3.044e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.9e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.553e-04
  |ΔD|_∞ = 6.088e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.9e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.107e-05
  |ΔD|_∞ = 1.218e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.9e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.021e-05
  |ΔD|_∞ = 2.435e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.9e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.043e-06
  |ΔD|_∞ = 4.871e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.4663e+00 J
  → Fracture energy : 4.3967e+00 J
  → Total energy    : 1.3863e+01 J


## Step 147/151: t = 3.50e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 146 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.92e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.161e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.344e-03
  |ΔD|_∞ = 1.478e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.92e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.269e-03
  |ΔD|_∞ = 2.956e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.92e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.538e-04
  |ΔD|_∞ = 5.913e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.92e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.075e-05
  |ΔD|_∞ = 1.183e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.92e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.015e-05
  |ΔD|_∞ = 2.365e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.92e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.030e-06
  |ΔD|_∞ = 4.730e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.4839e+00 J
  → Fracture energy : 4.4554e+00 J
  → Total energy    : 1.3939e+01 J


## Step 148/151: t = 3.53e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 147 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.94e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.066e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.257e-03
  |ΔD|_∞ = 1.511e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.94e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.251e-03
  |ΔD|_∞ = 3.022e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.94e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.503e-04
  |ΔD|_∞ = 6.045e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.94e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.005e-05
  |ΔD|_∞ = 1.209e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.94e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.001e-05
  |ΔD|_∞ = 2.418e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.94e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.002e-06
  |ΔD|_∞ = 4.836e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.5015e+00 J
  → Fracture energy : 4.5148e+00 J
  → Total energy    : 1.4016e+01 J


## Step 149/151: t = 3.55e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 148 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.96e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.947e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.158e-03
  |ΔD|_∞ = 1.537e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.96e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.232e-03
  |ΔD|_∞ = 3.074e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.96e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.463e-04
  |ΔD|_∞ = 6.148e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.96e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.973e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.926e-05
  |ΔD|_∞ = 1.230e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.96e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.852e-06
  |ΔD|_∞ = 2.459e-04

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.5204e+00 J
  → Fracture energy : 4.5745e+00 J
  → Total energy    : 1.4095e+01 J


## Step 150/151: t = 3.58e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 149 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.98e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.881e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.079e-03
  |ΔD|_∞ = 1.449e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.98e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.998e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.216e-03
  |ΔD|_∞ = 2.898e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.98e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.431e-04
  |ΔD|_∞ = 5.797e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.98e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.863e-05
  |ΔD|_∞ = 1.159e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [2.98e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.726e-06
  |ΔD|_∞ = 2.319e-04

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.5390e+00 J
  → Fracture energy : 4.6330e+00 J
  → Total energy    : 1.4172e+01 J


## Step 151/151: t = 3.60e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 150 | dt = 2.40e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.837e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.036e-03
  |ΔD|_∞ = 1.463e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.207e-03
  |ΔD|_∞ = 2.925e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.269e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.414e-04
  |ΔD|_∞ = 5.851e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.631e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.828e-05
  |ΔD|_∞ = 1.170e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [3e-05, 0.0]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.631e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.657e-06
  |ΔD|_∞ = 2.340e-04

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.5555e+00 J
  → Fracture energy : 4.6915e+00 J
  → Total energy    : 1.4247e+01 J

Simulation completed in 1550.09 s
Total time steps solved: 151
