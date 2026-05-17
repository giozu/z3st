Info    : Reading 'mesh.msh'...
Info    : 20594 nodes
Info    : 41516 elements
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
  → Time steps          : 200
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
  **[INFO]** Step-dependent Dirichlet list (2D), length 200
  **[INFO]** Dirichlet mechanical BC on 'steel' → [0.0, 4e-08] at region 'ymax'

Setting damage boundary conditions...
  **[INFO]** Dirichlet damage BC on 'steel' → D = 1.0 at region 'crack'


## Step 01/200: t = 0.00e+00 s | LHR = 0.00e+00 W/m

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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4e-08]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.000e+00
  |ΔD|_∞ = 8.000e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4e-08]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.794e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.000e-01
  |ΔD|_∞ = 1.600e-01

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4e-08]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.689e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.000e-02
  |ΔD|_∞ = 3.200e-02

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4e-08]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.751e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.000e-03
  |ΔD|_∞ = 6.400e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4e-08]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.845e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.600e-03
  |ΔD|_∞ = 1.280e-03

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4e-08]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.845e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.200e-04
  |ΔD|_∞ = 2.560e-04

Convergence check


#### Iteration 7/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4e-08]
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4e-08]
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4e-08]
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
  → Elastic energy  : 1.9304e-04 J
  → Fracture energy : 1.3669e+00 J
  → Total energy    : 1.3671e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0000.vtu


## Step 02/200: t = 1.81e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 1 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 8e-08]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.033e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.922e-05
  |ΔD|_∞ = 2.498e-05

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 8e-08]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.234e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.843e-06
  |ΔD|_∞ = 4.997e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.3793e-04 J
  → Fracture energy : 1.3669e+00 J
  → Total energy    : 1.3677e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0001.vtu


## Step 03/200: t = 3.62e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 2 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.2e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.333e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.531e-05
  |ΔD|_∞ = 3.206e-05

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.2e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.448e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.062e-06
  |ΔD|_∞ = 6.411e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6603e-03 J
  → Fracture energy : 1.3669e+00 J
  → Total energy    : 1.3686e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0002.vtu


## Step 04/200: t = 5.43e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 3 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.6e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.500e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.538e-05
  |ΔD|_∞ = 4.478e-05

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.6e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.277e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.076e-06
  |ΔD|_∞ = 8.955e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9516e-03 J
  → Fracture energy : 1.3669e+00 J
  → Total energy    : 1.3699e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0003.vtu


## Step 05/200: t = 7.24e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 4 | dt = 1.81e+01 s
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
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.000e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.562e-05
  |ΔD|_∞ = 5.774e-05

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.124e-06
  |ΔD|_∞ = 1.155e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.6118e-03 J
  → Fracture energy : 1.3669e+00 J
  → Total energy    : 1.3716e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0004.vtu


## Step 06/200: t = 9.05e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 5 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.4e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.667e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.588e-05
  |ΔD|_∞ = 7.075e-05

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.4e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.108e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.118e-05
  |ΔD|_∞ = 1.415e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.4e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.823e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.235e-06
  |ΔD|_∞ = 2.830e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.6408e-03 J
  → Fracture energy : 1.3669e+00 J
  → Total energy    : 1.3736e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0005.vtu


## Step 07/200: t = 1.09e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 6 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.8e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.429e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.437e-05
  |ΔD|_∞ = 8.153e-05

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.8e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.584e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.287e-05
  |ΔD|_∞ = 1.631e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.8e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.132e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.575e-06
  |ΔD|_∞ = 3.261e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.0386e-03 J
  → Fracture energy : 1.3669e+00 J
  → Total energy    : 1.3760e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0006.vtu


## Step 08/200: t = 1.27e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 7 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.2e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.250e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.432e-05
  |ΔD|_∞ = 9.418e-05

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.2e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.851e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.486e-05
  |ΔD|_∞ = 1.884e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.2e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.874e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.973e-06
  |ΔD|_∞ = 3.767e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1805e-02 J
  → Fracture energy : 1.3669e+00 J
  → Total energy    : 1.3787e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0007.vtu


## Step 09/200: t = 1.45e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 8 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.6e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.111e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.431e-05
  |ΔD|_∞ = 1.069e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.6e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.465e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.686e-05
  |ΔD|_∞ = 2.138e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.6e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.465e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.372e-06
  |ΔD|_∞ = 4.276e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4940e-02 J
  → Fracture energy : 1.3669e+00 J
  → Total energy    : 1.3819e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0008.vtu


## Step 10/200: t = 1.63e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 9 | dt = 1.81e+01 s
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
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.000e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.433e-05
  |ΔD|_∞ = 1.197e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.287e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.887e-05
  |ΔD|_∞ = 2.394e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.987e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.773e-06
  |ΔD|_∞ = 4.787e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8444e-02 J
  → Fracture energy : 1.3669e+00 J
  → Total energy    : 1.3854e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0009.vtu


## Step 11/200: t = 1.81e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 10 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.4e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.091e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.044e-04
  |ΔD|_∞ = 1.325e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.4e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.088e-05
  |ΔD|_∞ = 2.651e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.4e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.175e-06
  |ΔD|_∞ = 5.301e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2316e-02 J
  → Fracture energy : 1.3669e+00 J
  → Total energy    : 1.3893e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0010.vtu


## Step 12/200: t = 1.99e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 11 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.8e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.334e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.145e-04
  |ΔD|_∞ = 1.454e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.8e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.289e-05
  |ΔD|_∞ = 2.909e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.8e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.579e-06
  |ΔD|_∞ = 5.818e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6556e-02 J
  → Fracture energy : 1.3669e+00 J
  → Total energy    : 1.3935e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0011.vtu


## Step 13/200: t = 2.17e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 12 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.2e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.693e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.246e-04
  |ΔD|_∞ = 1.585e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.2e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.492e-05
  |ΔD|_∞ = 3.169e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.2e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.984e-06
  |ΔD|_∞ = 6.338e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1165e-02 J
  → Fracture energy : 1.3669e+00 J
  → Total energy    : 1.3981e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0012.vtu


## Step 14/200: t = 2.35e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 13 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.6e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.143e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.348e-04
  |ΔD|_∞ = 1.715e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.6e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.145e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.695e-05
  |ΔD|_∞ = 3.431e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.6e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.761e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.391e-06
  |ΔD|_∞ = 6.862e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.6141e-02 J
  → Fracture energy : 1.3670e+00 J
  → Total energy    : 1.4031e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0013.vtu


## Step 15/200: t = 2.53e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 14 | dt = 1.81e+01 s
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
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.667e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.450e-04
  |ΔD|_∞ = 1.847e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.900e-05
  |ΔD|_∞ = 3.695e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.039e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.800e-06
  |ΔD|_∞ = 7.390e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.1486e-02 J
  → Fracture energy : 1.3670e+00 J
  → Total energy    : 1.4084e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0014.vtu


## Step 16/200: t = 2.71e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 15 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.4e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.250e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.553e-04
  |ΔD|_∞ = 1.981e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.4e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.905e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.105e-05
  |ΔD|_∞ = 3.961e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.4e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.279e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.211e-06
  |ΔD|_∞ = 7.922e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.7198e-02 J
  → Fracture energy : 1.3670e+00 J
  → Total energy    : 1.4142e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0015.vtu


## Step 17/200: t = 2.89e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 16 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.8e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.883e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.656e-04
  |ΔD|_∞ = 2.115e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.8e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.990e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.312e-05
  |ΔD|_∞ = 4.230e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.8e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.624e-06
  |ΔD|_∞ = 8.459e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.3278e-02 J
  → Fracture energy : 1.3670e+00 J
  → Total energy    : 1.4202e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0016.vtu


## Step 18/200: t = 3.08e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 17 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.2e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.556e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.760e-04
  |ΔD|_∞ = 2.250e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.2e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.067e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.520e-05
  |ΔD|_∞ = 4.501e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.2e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.040e-06
  |ΔD|_∞ = 9.001e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.9725e-02 J
  → Fracture energy : 1.3670e+00 J
  → Total energy    : 1.4267e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0017.vtu


## Step 19/200: t = 3.26e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 18 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.6e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.264e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.865e-04
  |ΔD|_∞ = 2.387e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.6e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.729e-05
  |ΔD|_∞ = 4.775e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.6e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.459e-06
  |ΔD|_∞ = 9.549e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.6539e-02 J
  → Fracture energy : 1.3670e+00 J
  → Total energy    : 1.4335e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0018.vtu


## Step 20/200: t = 3.44e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 19 | dt = 1.81e+01 s
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
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.001e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.970e-04
  |ΔD|_∞ = 2.526e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 8e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.804e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.940e-05
  |ΔD|_∞ = 5.051e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 8e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.880e-06
  |ΔD|_∞ = 1.010e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.3720e-02 J
  → Fracture energy : 1.3670e+00 J
  → Total energy    : 1.4407e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0019.vtu


## Step 21/200: t = 3.62e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 20 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 8.4e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.762e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.076e-04
  |ΔD|_∞ = 2.666e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 8.4e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.759e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.152e-05
  |ΔD|_∞ = 5.332e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 8.4e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.759e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.305e-06
  |ΔD|_∞ = 1.066e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.1268e-02 J
  → Fracture energy : 1.3670e+00 J
  → Total energy    : 1.4483e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0020.vtu


## Step 22/200: t = 3.80e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 21 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 8.8e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.546e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.183e-04
  |ΔD|_∞ = 2.808e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 8.8e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.366e-05
  |ΔD|_∞ = 5.615e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 8.8e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.573e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.732e-06
  |ΔD|_∞ = 1.123e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.9183e-02 J
  → Fracture energy : 1.3670e+00 J
  → Total energy    : 1.4562e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0021.vtu


## Step 23/200: t = 3.98e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 22 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 9.2e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.348e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.291e-04
  |ΔD|_∞ = 2.951e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 9.2e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.582e-05
  |ΔD|_∞ = 5.902e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 9.2e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.164e-06
  |ΔD|_∞ = 1.180e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.7463e-02 J
  → Fracture energy : 1.3670e+00 J
  → Total energy    : 1.4645e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0022.vtu


## Step 24/200: t = 4.16e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 23 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 9.6e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.167e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.400e-04
  |ΔD|_∞ = 3.097e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 9.6e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.543e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.799e-05
  |ΔD|_∞ = 6.193e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 9.6e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.599e-06
  |ΔD|_∞ = 1.239e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0611e-01 J
  → Fracture energy : 1.3670e+00 J
  → Total energy    : 1.4731e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0023.vtu


## Step 25/200: t = 4.34e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 24 | dt = 1.81e+01 s
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
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.001e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.510e-04
  |ΔD|_∞ = 3.244e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.019e-05
  |ΔD|_∞ = 6.488e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.191e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.004e-05
  |ΔD|_∞ = 1.298e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.422e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.008e-06
  |ΔD|_∞ = 2.595e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1512e-01 J
  → Fracture energy : 1.3670e+00 J
  → Total energy    : 1.4822e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0024.vtu


## Step 26/200: t = 4.52e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 25 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.847e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.605e-04
  |ΔD|_∞ = 3.374e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.210e-05
  |ΔD|_∞ = 6.748e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.405e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.042e-05
  |ΔD|_∞ = 1.350e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.396e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.084e-06
  |ΔD|_∞ = 2.699e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2450e-01 J
  → Fracture energy : 1.3670e+00 J
  → Total energy    : 1.4915e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0025.vtu


## Step 27/200: t = 4.70e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 26 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.704e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.716e-04
  |ΔD|_∞ = 3.525e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.094e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.431e-05
  |ΔD|_∞ = 7.049e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.094e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.086e-05
  |ΔD|_∞ = 1.410e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.533e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.173e-06
  |ΔD|_∞ = 2.820e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3424e-01 J
  → Fracture energy : 1.3671e+00 J
  → Total energy    : 1.5013e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0026.vtu


## Step 28/200: t = 4.88e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 27 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.572e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.828e-04
  |ΔD|_∞ = 3.678e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.567e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.656e-05
  |ΔD|_∞ = 7.357e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.245e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.131e-05
  |ΔD|_∞ = 1.471e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.262e-06
  |ΔD|_∞ = 2.943e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4435e-01 J
  → Fracture energy : 1.3671e+00 J
  → Total energy    : 1.5114e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0027.vtu


## Step 29/200: t = 5.07e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 28 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.449e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.942e-04
  |ΔD|_∞ = 3.835e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.441e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.884e-05
  |ΔD|_∞ = 7.669e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.187e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.177e-05
  |ΔD|_∞ = 1.534e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.187e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.353e-06
  |ΔD|_∞ = 3.068e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5482e-01 J
  → Fracture energy : 1.3671e+00 J
  → Total energy    : 1.5219e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0028.vtu


## Step 30/200: t = 5.25e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 29 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.334e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.057e-04
  |ΔD|_∞ = 3.994e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.932e-20

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.114e-05
  |ΔD|_∞ = 7.987e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.842e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.223e-05
  |ΔD|_∞ = 1.597e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.842e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.446e-06
  |ΔD|_∞ = 3.195e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6566e-01 J
  → Fracture energy : 1.3671e+00 J
  → Total energy    : 1.5328e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0029.vtu


## Step 31/200: t = 5.43e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 30 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.227e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.173e-04
  |ΔD|_∞ = 4.156e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.375e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.347e-05
  |ΔD|_∞ = 8.312e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.224e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.269e-05
  |ΔD|_∞ = 1.662e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.418e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.539e-06
  |ΔD|_∞ = 3.325e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7686e-01 J
  → Fracture energy : 1.3672e+00 J
  → Total energy    : 1.5440e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0030.vtu


## Step 32/200: t = 5.61e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 31 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.126e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.291e-04
  |ΔD|_∞ = 4.321e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.207e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.583e-05
  |ΔD|_∞ = 8.642e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.016e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.317e-05
  |ΔD|_∞ = 1.728e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.633e-06
  |ΔD|_∞ = 3.457e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8842e-01 J
  → Fracture energy : 1.3672e+00 J
  → Total energy    : 1.5556e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0031.vtu


## Step 33/200: t = 5.79e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 32 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.031e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.411e-04
  |ΔD|_∞ = 4.490e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.675e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.822e-05
  |ΔD|_∞ = 8.979e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.096e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.364e-05
  |ΔD|_∞ = 1.796e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.729e-06
  |ΔD|_∞ = 3.592e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0034e-01 J
  → Fracture energy : 1.3672e+00 J
  → Total energy    : 1.5676e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0032.vtu


## Step 34/200: t = 5.97e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 33 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.942e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.532e-04
  |ΔD|_∞ = 4.662e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.273e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.065e-05
  |ΔD|_∞ = 9.324e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.526e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.413e-05
  |ΔD|_∞ = 1.865e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.592e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.826e-06
  |ΔD|_∞ = 3.730e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1263e-01 J
  → Fracture energy : 1.3673e+00 J
  → Total energy    : 1.5799e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0033.vtu


## Step 35/200: t = 6.15e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 34 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.858e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.656e-04
  |ΔD|_∞ = 4.838e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.311e-05
  |ΔD|_∞ = 9.675e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.301e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.462e-05
  |ΔD|_∞ = 1.935e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.301e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.925e-06
  |ΔD|_∞ = 3.870e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2528e-01 J
  → Fracture energy : 1.3673e+00 J
  → Total energy    : 1.5926e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0034.vtu


## Step 36/200: t = 6.33e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 35 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.779e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.781e-04
  |ΔD|_∞ = 5.018e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.939e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.561e-05
  |ΔD|_∞ = 1.004e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.303e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.512e-05
  |ΔD|_∞ = 2.007e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.212e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.024e-06
  |ΔD|_∞ = 4.014e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3829e-01 J
  → Fracture energy : 1.3673e+00 J
  → Total energy    : 1.6056e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0035.vtu


## Step 37/200: t = 6.51e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 36 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.704e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.908e-04
  |ΔD|_∞ = 5.202e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.815e-05
  |ΔD|_∞ = 1.040e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.563e-05
  |ΔD|_∞ = 2.081e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.126e-06
  |ΔD|_∞ = 4.162e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5167e-01 J
  → Fracture energy : 1.3674e+00 J
  → Total energy    : 1.6190e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0036.vtu


## Step 38/200: t = 6.69e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 37 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.633e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.037e-04
  |ΔD|_∞ = 5.391e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.295e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.074e-05
  |ΔD|_∞ = 1.078e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.295e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.615e-05
  |ΔD|_∞ = 2.156e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.229e-06
  |ΔD|_∞ = 4.313e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6540e-01 J
  → Fracture energy : 1.3674e+00 J
  → Total energy    : 1.6328e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0037.vtu


## Step 39/200: t = 6.87e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 38 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.565e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.168e-04
  |ΔD|_∞ = 5.584e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.186e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.337e-05
  |ΔD|_∞ = 1.117e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.667e-05
  |ΔD|_∞ = 2.234e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.702e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.335e-06
  |ΔD|_∞ = 4.467e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7949e-01 J
  → Fracture energy : 1.3675e+00 J
  → Total energy    : 1.6470e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0038.vtu


## Step 40/200: t = 7.06e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 39 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.501e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.302e-04
  |ΔD|_∞ = 5.783e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.079e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.605e-05
  |ΔD|_∞ = 1.157e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.521e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.721e-05
  |ΔD|_∞ = 2.313e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.079e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.442e-06
  |ΔD|_∞ = 4.626e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9394e-01 J
  → Fracture energy : 1.3675e+00 J
  → Total energy    : 1.6615e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0039.vtu


## Step 41/200: t = 7.24e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 40 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.440e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.439e-04
  |ΔD|_∞ = 5.987e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.961e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.878e-05
  |ΔD|_∞ = 1.197e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.233e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.776e-05
  |ΔD|_∞ = 2.395e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.551e-06
  |ΔD|_∞ = 4.790e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.0875e-01 J
  → Fracture energy : 1.3676e+00 J
  → Total energy    : 1.6764e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0040.vtu


## Step 42/200: t = 7.42e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 41 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.382e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.577e-04
  |ΔD|_∞ = 6.197e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.154e-05
  |ΔD|_∞ = 1.239e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.159e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.831e-05
  |ΔD|_∞ = 2.479e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.953e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.662e-06
  |ΔD|_∞ = 4.958e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.2392e-01 J
  → Fracture energy : 1.3677e+00 J
  → Total energy    : 1.6916e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0041.vtu


## Step 43/200: t = 7.60e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 42 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.327e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.720e-04
  |ΔD|_∞ = 6.413e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.467e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.440e-05
  |ΔD|_∞ = 1.283e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.888e-05
  |ΔD|_∞ = 2.565e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.568e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.776e-06
  |ΔD|_∞ = 5.131e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3945e-01 J
  → Fracture energy : 1.3678e+00 J
  → Total energy    : 1.7072e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0042.vtu


## Step 44/200: t = 7.78e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 43 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.274e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.865e-04
  |ΔD|_∞ = 6.636e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.444e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.731e-05
  |ΔD|_∞ = 1.327e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.223e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.946e-05
  |ΔD|_∞ = 2.654e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.223e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.892e-06
  |ΔD|_∞ = 5.309e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.5533e-01 J
  → Fracture energy : 1.3679e+00 J
  → Total energy    : 1.7232e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0043.vtu


## Step 45/200: t = 7.96e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 44 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.224e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.014e-04
  |ΔD|_∞ = 6.866e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.837e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.003e-04
  |ΔD|_∞ = 1.373e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.963e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.006e-05
  |ΔD|_∞ = 2.746e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.011e-06
  |ΔD|_∞ = 5.493e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.7157e-01 J
  → Fracture energy : 1.3679e+00 J
  → Total energy    : 1.7395e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0044.vtu


## Step 46/200: t = 8.14e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 45 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.176e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.166e-04
  |ΔD|_∞ = 7.103e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.348e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.033e-04
  |ΔD|_∞ = 1.421e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.348e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.066e-05
  |ΔD|_∞ = 2.841e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.348e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.133e-06
  |ΔD|_∞ = 5.683e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.8816e-01 J
  → Fracture energy : 1.3680e+00 J
  → Total energy    : 1.7562e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0045.vtu


## Step 47/200: t = 8.32e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 46 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.130e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.322e-04
  |ΔD|_∞ = 7.349e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.608e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.064e-04
  |ΔD|_∞ = 1.470e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.609e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.129e-05
  |ΔD|_∞ = 2.940e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.386e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.257e-06
  |ΔD|_∞ = 5.879e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.0511e-01 J
  → Fracture energy : 1.3681e+00 J
  → Total energy    : 1.7733e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0046.vtu


## Step 48/200: t = 8.50e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 47 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.085e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.482e-04
  |ΔD|_∞ = 7.603e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.110e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.096e-04
  |ΔD|_∞ = 1.521e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.912e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.193e-05
  |ΔD|_∞ = 3.041e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.195e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.385e-06
  |ΔD|_∞ = 6.083e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.2241e-01 J
  → Fracture energy : 1.3683e+00 J
  → Total energy    : 1.7907e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0047.vtu


## Step 49/200: t = 8.68e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 48 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.043e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.646e-04
  |ΔD|_∞ = 7.867e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.178e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.129e-04
  |ΔD|_∞ = 1.573e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.366e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.258e-05
  |ΔD|_∞ = 3.147e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.923e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.517e-06
  |ΔD|_∞ = 6.293e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.4006e-01 J
  → Fracture energy : 1.3684e+00 J
  → Total energy    : 1.8084e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0048.vtu


## Step 50/200: t = 8.86e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 49 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.002e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.815e-04
  |ΔD|_∞ = 8.140e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.475e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.163e-04
  |ΔD|_∞ = 1.628e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.230e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.326e-05
  |ΔD|_∞ = 3.256e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.177e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.652e-06
  |ΔD|_∞ = 6.512e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.5807e-01 J
  → Fracture energy : 1.3685e+00 J
  → Total energy    : 1.8266e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0049.vtu


## Step 51/200: t = 9.05e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 50 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.963e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.988e-04
  |ΔD|_∞ = 8.425e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.392e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.198e-04
  |ΔD|_∞ = 1.685e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.048e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.395e-05
  |ΔD|_∞ = 3.370e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.649e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.791e-06
  |ΔD|_∞ = 6.740e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.7642e-01 J
  → Fracture energy : 1.3687e+00 J
  → Total energy    : 1.8451e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0050.vtu


## Step 52/200: t = 9.23e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 51 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.925e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.168e-04
  |ΔD|_∞ = 8.721e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.529e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.234e-04
  |ΔD|_∞ = 1.744e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.530e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.467e-05
  |ΔD|_∞ = 3.488e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.070e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.934e-06
  |ΔD|_∞ = 6.977e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.9513e-01 J
  → Fracture energy : 1.3688e+00 J
  → Total energy    : 1.8639e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0051.vtu


## Step 53/200: t = 9.41e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 52 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.889e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.352e-04
  |ΔD|_∞ = 9.029e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.316e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.270e-04
  |ΔD|_∞ = 1.806e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.764e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.541e-05
  |ΔD|_∞ = 3.612e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.316e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.082e-06
  |ΔD|_∞ = 7.223e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.1418e-01 J
  → Fracture energy : 1.3690e+00 J
  → Total energy    : 1.8832e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0052.vtu


## Step 54/200: t = 9.59e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 53 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.854e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.544e-04
  |ΔD|_∞ = 9.352e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.699e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.309e-04
  |ΔD|_∞ = 1.870e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.789e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.617e-05
  |ΔD|_∞ = 3.741e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.114e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.235e-06
  |ΔD|_∞ = 7.481e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.3358e-01 J
  → Fracture energy : 1.3691e+00 J
  → Total energy    : 1.9027e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0053.vtu


## Step 55/200: t = 9.77e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 54 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.821e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.741e-04
  |ΔD|_∞ = 9.689e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.348e-04
  |ΔD|_∞ = 1.938e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.697e-05
  |ΔD|_∞ = 3.875e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.081e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.393e-06
  |ΔD|_∞ = 7.751e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.5333e-01 J
  → Fracture energy : 1.3693e+00 J
  → Total energy    : 1.9227e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0054.vtu


## Step 56/200: t = 9.95e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 55 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.789e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.947e-04
  |ΔD|_∞ = 1.004e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.389e-04
  |ΔD|_∞ = 2.008e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.336e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.779e-05
  |ΔD|_∞ = 4.017e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.858e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.557e-06
  |ΔD|_∞ = 8.033e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.7342e-01 J
  → Fracture energy : 1.3695e+00 J
  → Total energy    : 1.9430e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0055.vtu


## Step 57/200: t = 1.01e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 56 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.757e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.160e-04
  |ΔD|_∞ = 1.041e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.538e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.432e-04
  |ΔD|_∞ = 2.082e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.864e-05
  |ΔD|_∞ = 4.165e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.538e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.728e-06
  |ΔD|_∞ = 8.330e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.9386e-01 J
  → Fracture energy : 1.3698e+00 J
  → Total energy    : 1.9636e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0056.vtu


## Step 58/200: t = 1.03e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 57 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.727e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.381e-04
  |ΔD|_∞ = 1.080e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.799e-20

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.476e-04
  |ΔD|_∞ = 2.160e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.371e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.953e-05
  |ΔD|_∞ = 4.321e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.370e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.905e-06
  |ΔD|_∞ = 8.642e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.1463e-01 J
  → Fracture energy : 1.3700e+00 J
  → Total energy    : 1.9846e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0057.vtu


## Step 59/200: t = 1.05e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 58 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.698e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.612e-04
  |ΔD|_∞ = 1.121e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.715e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.522e-04
  |ΔD|_∞ = 2.243e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.410e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.045e-05
  |ΔD|_∞ = 4.485e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.502e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.090e-06
  |ΔD|_∞ = 8.970e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.3575e-01 J
  → Fracture energy : 1.3702e+00 J
  → Total energy    : 2.0060e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0058.vtu


## Step 60/200: t = 1.07e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 59 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.670e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.854e-04
  |ΔD|_∞ = 1.166e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.571e-04
  |ΔD|_∞ = 2.333e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.739e-20

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.142e-05
  |ΔD|_∞ = 4.665e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.283e-06
  |ΔD|_∞ = 9.330e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.5720e-01 J
  → Fracture energy : 1.3705e+00 J
  → Total energy    : 2.0277e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0059.vtu


## Step 61/200: t = 1.09e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 60 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.643e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.107e-04
  |ΔD|_∞ = 1.215e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.538e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.621e-04
  |ΔD|_∞ = 2.430e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.393e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.243e-05
  |ΔD|_∞ = 4.860e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.559e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.485e-06
  |ΔD|_∞ = 9.720e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.7899e-01 J
  → Fracture energy : 1.3708e+00 J
  → Total energy    : 2.0498e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0060.vtu


## Step 62/200: t = 1.10e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 61 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.617e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.372e-04
  |ΔD|_∞ = 1.267e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.374e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.674e-04
  |ΔD|_∞ = 2.534e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.374e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.349e-05
  |ΔD|_∞ = 5.068e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.550e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.698e-06
  |ΔD|_∞ = 1.014e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.0111e-01 J
  → Fracture energy : 1.3711e+00 J
  → Total energy    : 2.0722e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0061.vtu


## Step 63/200: t = 1.12e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 62 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.591e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.652e-04
  |ΔD|_∞ = 1.322e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.699e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.730e-04
  |ΔD|_∞ = 2.645e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.003e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.461e-05
  |ΔD|_∞ = 5.290e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.246e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.921e-06
  |ΔD|_∞ = 1.058e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.2356e-01 J
  → Fracture energy : 1.3714e+00 J
  → Total energy    : 2.0950e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0062.vtu


## Step 64/200: t = 1.14e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 63 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.567e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.946e-04
  |ΔD|_∞ = 1.382e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.424e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.789e-04
  |ΔD|_∞ = 2.764e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.426e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.578e-05
  |ΔD|_∞ = 5.529e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.157e-06
  |ΔD|_∞ = 1.106e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.4635e-01 J
  → Fracture energy : 1.3718e+00 J
  → Total energy    : 2.1181e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0063.vtu


## Step 65/200: t = 1.16e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 64 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.543e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.258e-04
  |ΔD|_∞ = 1.447e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.234e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.852e-04
  |ΔD|_∞ = 2.893e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.661e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.703e-05
  |ΔD|_∞ = 5.786e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.722e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.406e-06
  |ΔD|_∞ = 1.157e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.6946e-01 J
  → Fracture energy : 1.3722e+00 J
  → Total energy    : 2.1416e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0064.vtu


## Step 66/200: t = 1.18e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 65 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.520e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.589e-04
  |ΔD|_∞ = 1.516e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.458e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.918e-04
  |ΔD|_∞ = 3.032e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.459e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.835e-05
  |ΔD|_∞ = 6.065e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.418e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.671e-06
  |ΔD|_∞ = 1.213e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.9289e-01 J
  → Fracture energy : 1.3726e+00 J
  → Total energy    : 2.1654e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0065.vtu


## Step 67/200: t = 1.19e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 66 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.498e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.942e-04
  |ΔD|_∞ = 1.592e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.592e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.988e-04
  |ΔD|_∞ = 3.184e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.337e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.977e-05
  |ΔD|_∞ = 6.367e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.953e-06
  |ΔD|_∞ = 1.273e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.1664e-01 J
  → Fracture energy : 1.3730e+00 J
  → Total energy    : 2.1896e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0066.vtu


## Step 68/200: t = 1.21e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 67 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.476e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.032e-03
  |ΔD|_∞ = 1.674e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.143e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.064e-04
  |ΔD|_∞ = 3.349e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.143e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.128e-05
  |ΔD|_∞ = 6.698e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.726e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.255e-06
  |ΔD|_∞ = 1.340e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.4070e-01 J
  → Fracture energy : 1.3735e+00 J
  → Total energy    : 2.2142e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0067.vtu


## Step 69/200: t = 1.23e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 68 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.455e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.072e-03
  |ΔD|_∞ = 1.765e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.461e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.145e-04
  |ΔD|_∞ = 3.530e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.290e-05
  |ΔD|_∞ = 7.060e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.580e-06
  |ΔD|_∞ = 1.412e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.6508e-01 J
  → Fracture energy : 1.3740e+00 J
  → Total energy    : 2.2391e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0068.vtu


## Step 70/200: t = 1.25e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 69 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.435e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.116e-03
  |ΔD|_∞ = 1.865e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.793e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.233e-04
  |ΔD|_∞ = 3.730e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.465e-05
  |ΔD|_∞ = 7.460e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.793e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.930e-06
  |ΔD|_∞ = 1.492e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.8976e-01 J
  → Fracture energy : 1.3745e+00 J
  → Total energy    : 2.2643e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0069.vtu


## Step 71/200: t = 1.27e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 70 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.415e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.164e-03
  |ΔD|_∞ = 1.976e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.877e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.328e-04
  |ΔD|_∞ = 3.952e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.877e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.655e-05
  |ΔD|_∞ = 7.903e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.326e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.311e-06
  |ΔD|_∞ = 1.581e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.1474e-01 J
  → Fracture energy : 1.3751e+00 J
  → Total energy    : 2.2899e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0070.vtu


## Step 72/200: t = 1.28e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 71 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.396e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.216e-03
  |ΔD|_∞ = 2.100e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.302e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.432e-04
  |ΔD|_∞ = 4.199e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.678e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.863e-05
  |ΔD|_∞ = 8.399e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.679e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.726e-06
  |ΔD|_∞ = 1.680e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.4002e-01 J
  → Fracture energy : 1.3758e+00 J
  → Total energy    : 2.3158e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0071.vtu


## Step 73/200: t = 1.30e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 72 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.378e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.273e-03
  |ΔD|_∞ = 2.239e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.444e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.546e-04
  |ΔD|_∞ = 4.478e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.354e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.092e-05
  |ΔD|_∞ = 8.956e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.355e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.018e-05
  |ΔD|_∞ = 1.791e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.885e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.037e-06
  |ΔD|_∞ = 3.583e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.6558e-01 J
  → Fracture energy : 1.3765e+00 J
  → Total energy    : 2.3420e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0072.vtu


## Step 74/200: t = 1.32e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 73 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.360e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.335e-03
  |ΔD|_∞ = 2.396e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.635e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.671e-04
  |ΔD|_∞ = 4.792e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.500e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.341e-05
  |ΔD|_∞ = 9.584e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.416e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.068e-05
  |ΔD|_∞ = 1.917e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.418e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.136e-06
  |ΔD|_∞ = 3.834e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.9141e-01 J
  → Fracture energy : 1.3772e+00 J
  → Total energy    : 2.3686e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0073.vtu


## Step 75/200: t = 1.34e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 74 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.343e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.406e-03
  |ΔD|_∞ = 2.577e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.849e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.811e-04
  |ΔD|_∞ = 5.154e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.849e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.623e-05
  |ΔD|_∞ = 1.031e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.109e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.125e-05
  |ΔD|_∞ = 2.061e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.109e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.249e-06
  |ΔD|_∞ = 4.123e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0175e+00 J
  → Fracture energy : 1.3780e+00 J
  → Total energy    : 2.3956e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0074.vtu


## Step 76/200: t = 1.36e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 75 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.327e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.485e-03
  |ΔD|_∞ = 2.786e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.694e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.971e-04
  |ΔD|_∞ = 5.573e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.208e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.942e-05
  |ΔD|_∞ = 1.115e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.294e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.188e-05
  |ΔD|_∞ = 2.229e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.377e-06
  |ΔD|_∞ = 4.458e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0439e+00 J
  → Fracture energy : 1.3790e+00 J
  → Total energy    : 2.4228e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0075.vtu


## Step 77/200: t = 1.37e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 76 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.311e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.577e-03
  |ΔD|_∞ = 3.033e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.221e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.153e-04
  |ΔD|_∞ = 6.065e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.475e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.307e-05
  |ΔD|_∞ = 1.213e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.475e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.261e-05
  |ΔD|_∞ = 2.426e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.523e-06
  |ΔD|_∞ = 4.852e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0705e+00 J
  → Fracture energy : 1.3800e+00 J
  → Total energy    : 2.4504e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0076.vtu


## Step 78/200: t = 1.39e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 77 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.296e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.682e-03
  |ΔD|_∞ = 3.325e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.171e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.365e-04
  |ΔD|_∞ = 6.651e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.399e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.730e-05
  |ΔD|_∞ = 1.330e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.461e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.346e-05
  |ΔD|_∞ = 2.660e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.461e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.692e-06
  |ΔD|_∞ = 5.320e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0973e+00 J
  → Fracture energy : 1.3811e+00 J
  → Total energy    : 2.4784e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0077.vtu


## Step 79/200: t = 1.41e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 78 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.281e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.807e-03
  |ΔD|_∞ = 3.690e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.699e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.614e-04
  |ΔD|_∞ = 7.380e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.370e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.228e-05
  |ΔD|_∞ = 1.476e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.370e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.446e-05
  |ΔD|_∞ = 2.952e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.685e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.891e-06
  |ΔD|_∞ = 5.904e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1243e+00 J
  → Fracture energy : 1.3823e+00 J
  → Total energy    : 2.5066e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0078.vtu


## Step 80/200: t = 1.43e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 79 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.268e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.956e-03
  |ΔD|_∞ = 4.156e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.358e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.913e-04
  |ΔD|_∞ = 8.311e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.083e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.826e-05
  |ΔD|_∞ = 1.662e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.734e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.565e-05
  |ΔD|_∞ = 3.325e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.073e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.130e-06
  |ΔD|_∞ = 6.649e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1514e+00 J
  → Fracture energy : 1.3837e+00 J
  → Total energy    : 2.5352e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0079.vtu


## Step 81/200: t = 1.45e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 80 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.256e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.140e-03
  |ΔD|_∞ = 4.741e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.280e-04
  |ΔD|_∞ = 9.481e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.837e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.559e-05
  |ΔD|_∞ = 1.896e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.919e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.712e-05
  |ΔD|_∞ = 3.792e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.777e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.424e-06
  |ΔD|_∞ = 7.585e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1787e+00 J
  → Fracture energy : 1.3853e+00 J
  → Total energy    : 2.5640e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0080.vtu


## Step 82/200: t = 1.47e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 81 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.245e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.371e-03
  |ΔD|_∞ = 5.491e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.048e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.742e-04
  |ΔD|_∞ = 1.098e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.241e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.483e-05
  |ΔD|_∞ = 2.197e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.722e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.897e-05
  |ΔD|_∞ = 4.393e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.148e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.793e-06
  |ΔD|_∞ = 8.786e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2060e+00 J
  → Fracture energy : 1.3872e+00 J
  → Total energy    : 2.5932e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0081.vtu


## Step 83/200: t = 1.48e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 82 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.237e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.670e-03
  |ΔD|_∞ = 6.475e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.768e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.341e-04
  |ΔD|_∞ = 1.295e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.768e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.068e-04
  |ΔD|_∞ = 2.590e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.136e-05
  |ΔD|_∞ = 5.180e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.570e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.272e-06
  |ΔD|_∞ = 1.036e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2332e+00 J
  → Fracture energy : 1.3894e+00 J
  → Total energy    : 2.6226e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0082.vtu


## Step 84/200: t = 1.50e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 83 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.231e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.072e-03
  |ΔD|_∞ = 7.793e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.491e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.143e-04
  |ΔD|_∞ = 1.559e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.490e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.229e-04
  |ΔD|_∞ = 3.117e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.457e-05
  |ΔD|_∞ = 6.235e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.490e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.915e-06
  |ΔD|_∞ = 1.247e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2602e+00 J
  → Fracture energy : 1.3921e+00 J
  → Total energy    : 2.6523e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0083.vtu


## Step 85/200: t = 1.52e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 84 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.232e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.631e-03
  |ΔD|_∞ = 9.578e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.877e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.261e-04
  |ΔD|_∞ = 1.916e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.450e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.452e-04
  |ΔD|_∞ = 3.831e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.904e-05
  |ΔD|_∞ = 7.663e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.809e-06
  |ΔD|_∞ = 1.533e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2867e+00 J
  → Fracture energy : 1.3954e+00 J
  → Total energy    : 2.6821e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0084.vtu


## Step 86/200: t = 1.54e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 85 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.244e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.439e-03
  |ΔD|_∞ = 1.230e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.747e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.879e-04
  |ΔD|_∞ = 2.459e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.776e-04
  |ΔD|_∞ = 4.918e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.747e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.552e-05
  |ΔD|_∞ = 9.836e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.782e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.103e-06
  |ΔD|_∞ = 1.967e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3123e+00 J
  → Fracture energy : 1.3998e+00 J
  → Total energy    : 2.7121e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0085.vtu


## Step 87/200: t = 1.56e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 86 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.279e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.649e-03
  |ΔD|_∞ = 1.666e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.012e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.130e-03
  |ΔD|_∞ = 3.333e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.009e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.259e-04
  |ΔD|_∞ = 6.665e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.519e-05
  |ΔD|_∞ = 1.333e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.009e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.038e-06
  |ΔD|_∞ = 2.666e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3360e+00 J
  → Fracture energy : 1.4058e+00 J
  → Total energy    : 2.7419e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0086.vtu


## Step 88/200: t = 1.57e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 87 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.371e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.489e-03
  |ΔD|_∞ = 2.360e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.498e-03
  |ΔD|_∞ = 4.720e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.070e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.996e-04
  |ΔD|_∞ = 9.441e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.072e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.991e-05
  |ΔD|_∞ = 1.888e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.400e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.198e-05
  |ΔD|_∞ = 3.776e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.396e-06
  |ΔD|_∞ = 7.552e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3562e+00 J
  → Fracture energy : 1.4147e+00 J
  → Total energy    : 2.7709e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0087.vtu


## Step 89/200: t = 1.59e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 88 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.612e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.024e-02
  |ΔD|_∞ = 3.497e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.048e-03
  |ΔD|_∞ = 6.993e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.096e-04
  |ΔD|_∞ = 1.399e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.196e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.192e-05
  |ΔD|_∞ = 2.797e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.628e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.638e-05
  |ΔD|_∞ = 5.595e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.628e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.277e-06
  |ΔD|_∞ = 1.119e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3690e+00 J
  → Fracture energy : 1.4287e+00 J
  → Total energy    : 2.7977e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0088.vtu


## Step 90/200: t = 1.61e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 89 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.248e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.397e-02
  |ΔD|_∞ = 5.270e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.794e-03
  |ΔD|_∞ = 1.054e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.569e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.587e-04
  |ΔD|_∞ = 2.108e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.636e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.117e-04
  |ΔD|_∞ = 4.216e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.302e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.235e-05
  |ΔD|_∞ = 8.432e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.302e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.470e-06
  |ΔD|_∞ = 1.686e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3667e+00 J
  → Fracture energy : 1.4515e+00 J
  → Total energy    : 2.8181e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0089.vtu


## Step 91/200: t = 1.63e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 90 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.763e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.796e-02
  |ΔD|_∞ = 6.054e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.141e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.592e-03
  |ΔD|_∞ = 1.211e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.173e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.183e-04
  |ΔD|_∞ = 2.422e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.848e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.437e-04
  |ΔD|_∞ = 4.843e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.873e-05
  |ΔD|_∞ = 9.687e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.846e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.747e-06
  |ΔD|_∞ = 1.937e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3362e+00 J
  → Fracture energy : 1.4864e+00 J
  → Total energy    : 2.8226e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0090.vtu


## Step 92/200: t = 1.65e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 91 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.282e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.034e-02
  |ΔD|_∞ = 6.521e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.132e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.067e-03
  |ΔD|_∞ = 1.304e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.839e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.135e-04
  |ΔD|_∞ = 2.608e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.837e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.627e-04
  |ΔD|_∞ = 5.217e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.140e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.254e-05
  |ΔD|_∞ = 1.043e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.703e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.508e-06
  |ΔD|_∞ = 2.087e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2709e+00 J
  → Fracture energy : 1.5305e+00 J
  → Total energy    : 2.8015e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0091.vtu


## Step 93/200: t = 1.66e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 92 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.591e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.999e-02
  |ΔD|_∞ = 8.428e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.098e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.997e-03
  |ΔD|_∞ = 1.686e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.098e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.995e-04
  |ΔD|_∞ = 3.371e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.531e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.599e-04
  |ΔD|_∞ = 6.742e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.021e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.198e-05
  |ΔD|_∞ = 1.348e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.021e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.396e-06
  |ΔD|_∞ = 2.697e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1868e+00 J
  → Fracture energy : 1.5731e+00 J
  → Total energy    : 2.7598e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0092.vtu


## Step 94/200: t = 1.68e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 93 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.789e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.769e-02
  |ΔD|_∞ = 1.063e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.267e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.537e-03
  |ΔD|_∞ = 2.127e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.266e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.075e-04
  |ΔD|_∞ = 4.254e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.268e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.415e-04
  |ΔD|_∞ = 8.507e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.268e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.830e-05
  |ΔD|_∞ = 1.701e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.660e-06
  |ΔD|_∞ = 3.403e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1174e+00 J
  → Fracture energy : 1.6072e+00 J
  → Total energy    : 2.7246e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0093.vtu


## Step 95/200: t = 1.70e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 94 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.332e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.386e-02
  |ΔD|_∞ = 1.469e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.950e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.771e-03
  |ΔD|_∞ = 2.939e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.738e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.543e-04
  |ΔD|_∞ = 5.878e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.670e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.109e-04
  |ΔD|_∞ = 1.176e-03

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.892e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.217e-05
  |ΔD|_∞ = 2.351e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.892e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.434e-06
  |ΔD|_∞ = 4.702e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0684e+00 J
  → Fracture energy : 1.6314e+00 J
  → Total energy    : 2.6998e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0094.vtu


## Step 96/200: t = 1.72e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 95 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.219e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.309e-03
  |ΔD|_∞ = 1.150e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.345e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.462e-03
  |ΔD|_∞ = 2.300e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.356e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.924e-04
  |ΔD|_∞ = 4.601e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.805e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.847e-05
  |ΔD|_∞ = 9.202e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.245e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.169e-05
  |ΔD|_∞ = 1.840e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.347e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.339e-06
  |ΔD|_∞ = 3.681e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0505e+00 J
  → Fracture energy : 1.6436e+00 J
  → Total energy    : 2.6940e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0095.vtu


## Step 97/200: t = 1.74e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 96 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.882e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.569e-03
  |ΔD|_∞ = 4.250e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.139e-04
  |ΔD|_∞ = 8.499e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.467e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.428e-04
  |ΔD|_∞ = 1.700e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.490e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.855e-05
  |ΔD|_∞ = 3.400e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.112e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.711e-06
  |ΔD|_∞ = 6.800e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0630e+00 J
  → Fracture energy : 1.6481e+00 J
  → Total energy    : 2.7112e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0096.vtu


## Step 98/200: t = 1.75e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 97 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.329e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.423e-03
  |ΔD|_∞ = 2.775e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.720e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.846e-04
  |ΔD|_∞ = 5.551e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.440e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.692e-05
  |ΔD|_∞ = 1.110e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.440e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.938e-05
  |ΔD|_∞ = 2.220e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.390e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.877e-06
  |ΔD|_∞ = 4.441e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0812e+00 J
  → Fracture energy : 1.6509e+00 J
  → Total energy    : 2.7321e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0097.vtu


## Step 99/200: t = 1.77e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 98 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.138e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.630e-03
  |ΔD|_∞ = 1.818e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.259e-04
  |ΔD|_∞ = 3.636e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.519e-05
  |ΔD|_∞ = 7.271e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.210e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.304e-05
  |ΔD|_∞ = 1.454e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.210e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.607e-06
  |ΔD|_∞ = 2.909e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1012e+00 J
  → Fracture energy : 1.6527e+00 J
  → Total energy    : 2.7539e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0098.vtu


## Step 100/200: t = 1.79e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 99 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.049e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.241e-03
  |ΔD|_∞ = 1.651e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.881e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.482e-04
  |ΔD|_∞ = 3.302e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.541e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.965e-05
  |ΔD|_∞ = 6.604e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.929e-06
  |ΔD|_∞ = 1.321e-04

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1221e+00 J
  → Fracture energy : 1.6541e+00 J
  → Total energy    : 2.7762e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0099.vtu


## Step 101/200: t = 1.81e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 100 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.019e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.073e-03
  |ΔD|_∞ = 1.503e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.285e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.147e-04
  |ΔD|_∞ = 3.006e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.285e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.293e-05
  |ΔD|_∞ = 6.013e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.265e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.586e-06
  |ΔD|_∞ = 1.203e-04

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1434e+00 J
  → Fracture energy : 1.6552e+00 J
  → Total energy    : 2.7986e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0100.vtu


## Step 102/200: t = 1.83e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 101 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.003e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.955e-04
  |ΔD|_∞ = 1.377e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.917e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.991e-04
  |ΔD|_∞ = 2.755e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.982e-05
  |ΔD|_∞ = 5.509e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.964e-06
  |ΔD|_∞ = 1.102e-04

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1650e+00 J
  → Fracture energy : 1.6563e+00 J
  → Total energy    : 2.8213e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0101.vtu


## Step 103/200: t = 1.85e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 102 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.905e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.546e-04
  |ΔD|_∞ = 1.266e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.909e-04
  |ΔD|_∞ = 2.532e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.220e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.818e-05
  |ΔD|_∞ = 5.065e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.466e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.637e-06
  |ΔD|_∞ = 1.013e-04

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1868e+00 J
  → Fracture energy : 1.6574e+00 J
  → Total energy    : 2.8442e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0102.vtu


## Step 104/200: t = 1.86e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 103 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.793e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.281e-04
  |ΔD|_∞ = 1.158e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.256e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.856e-04
  |ΔD|_∞ = 2.316e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.256e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.713e-05
  |ΔD|_∞ = 4.632e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.425e-06
  |ΔD|_∞ = 9.264e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2089e+00 J
  → Fracture energy : 1.6584e+00 J
  → Total energy    : 2.8673e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0103.vtu


## Step 105/200: t = 1.88e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 104 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.690e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.064e-04
  |ΔD|_∞ = 1.067e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.813e-04
  |ΔD|_∞ = 2.133e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.427e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.626e-05
  |ΔD|_∞ = 4.266e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.252e-06
  |ΔD|_∞ = 8.532e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2311e+00 J
  → Fracture energy : 1.6595e+00 J
  → Total energy    : 2.8906e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0104.vtu


## Step 106/200: t = 1.90e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 105 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.593e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.852e-04
  |ΔD|_∞ = 1.057e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.038e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.770e-04
  |ΔD|_∞ = 2.113e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.353e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.541e-05
  |ΔD|_∞ = 4.227e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.062e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.082e-06
  |ΔD|_∞ = 8.453e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2536e+00 J
  → Fracture energy : 1.6605e+00 J
  → Total energy    : 2.9141e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0105.vtu


## Step 107/200: t = 1.92e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 106 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.498e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.642e-04
  |ΔD|_∞ = 1.035e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.458e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.728e-04
  |ΔD|_∞ = 2.071e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.460e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.457e-05
  |ΔD|_∞ = 4.142e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.389e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.913e-06
  |ΔD|_∞ = 8.284e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2763e+00 J
  → Fracture energy : 1.6615e+00 J
  → Total energy    : 2.9379e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0106.vtu


## Step 108/200: t = 1.94e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 107 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.404e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.444e-04
  |ΔD|_∞ = 1.004e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.240e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.689e-04
  |ΔD|_∞ = 2.007e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.378e-05
  |ΔD|_∞ = 4.014e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.755e-06
  |ΔD|_∞ = 8.029e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2992e+00 J
  → Fracture energy : 1.6625e+00 J
  → Total energy    : 2.9618e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0107.vtu


## Step 109/200: t = 1.95e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 108 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.312e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.270e-04
  |ΔD|_∞ = 9.636e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.363e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.654e-04
  |ΔD|_∞ = 1.927e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.363e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.308e-05
  |ΔD|_∞ = 3.854e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.325e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.616e-06
  |ΔD|_∞ = 7.709e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3224e+00 J
  → Fracture energy : 1.6635e+00 J
  → Total energy    : 2.9859e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0108.vtu


## Step 110/200: t = 1.97e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 109 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.222e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.125e-04
  |ΔD|_∞ = 9.191e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.124e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.625e-04
  |ΔD|_∞ = 1.838e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.392e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.250e-05
  |ΔD|_∞ = 3.676e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.500e-06
  |ΔD|_∞ = 7.353e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3457e+00 J
  → Fracture energy : 1.6645e+00 J
  → Total energy    : 3.0102e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0109.vtu


## Step 111/200: t = 1.99e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 110 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.134e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.012e-04
  |ΔD|_∞ = 8.730e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.332e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.602e-04
  |ΔD|_∞ = 1.746e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.324e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.205e-05
  |ΔD|_∞ = 3.492e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.388e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.409e-06
  |ΔD|_∞ = 6.984e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3693e+00 J
  → Fracture energy : 1.6655e+00 J
  → Total energy    : 3.0348e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0110.vtu


## Step 112/200: t = 2.01e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 111 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.048e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.930e-04
  |ΔD|_∞ = 8.276e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.128e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.586e-04
  |ΔD|_∞ = 1.655e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.132e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.172e-05
  |ΔD|_∞ = 3.310e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.008e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.344e-06
  |ΔD|_∞ = 6.621e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3931e+00 J
  → Fracture energy : 1.6664e+00 J
  → Total energy    : 3.0595e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0111.vtu


## Step 113/200: t = 2.03e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 112 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.963e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.880e-04
  |ΔD|_∞ = 7.842e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.613e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.576e-04
  |ΔD|_∞ = 1.568e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.523e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.152e-05
  |ΔD|_∞ = 3.137e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.304e-06
  |ΔD|_∞ = 6.274e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4171e+00 J
  → Fracture energy : 1.6674e+00 J
  → Total energy    : 3.0845e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0112.vtu


## Step 114/200: t = 2.04e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 113 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.881e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.862e-04
  |ΔD|_∞ = 7.836e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.258e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.572e-04
  |ΔD|_∞ = 1.567e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.011e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.145e-05
  |ΔD|_∞ = 3.134e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.729e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.289e-06
  |ΔD|_∞ = 6.269e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4413e+00 J
  → Fracture energy : 1.6684e+00 J
  → Total energy    : 3.1096e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0113.vtu


## Step 115/200: t = 2.06e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 114 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.800e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.875e-04
  |ΔD|_∞ = 7.943e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.575e-04
  |ΔD|_∞ = 1.589e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.147e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.150e-05
  |ΔD|_∞ = 3.177e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.147e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.300e-06
  |ΔD|_∞ = 6.354e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4657e+00 J
  → Fracture energy : 1.6693e+00 J
  → Total energy    : 3.1350e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0114.vtu


## Step 116/200: t = 2.08e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 115 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.722e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.918e-04
  |ΔD|_∞ = 8.069e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.019e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.584e-04
  |ΔD|_∞ = 1.614e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.758e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.167e-05
  |ΔD|_∞ = 3.227e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.031e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.335e-06
  |ΔD|_∞ = 6.455e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4903e+00 J
  → Fracture energy : 1.6703e+00 J
  → Total energy    : 3.1605e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0115.vtu


## Step 117/200: t = 2.10e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 116 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.645e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.991e-04
  |ΔD|_∞ = 8.218e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.003e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.598e-04
  |ΔD|_∞ = 1.644e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.938e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.196e-05
  |ΔD|_∞ = 3.287e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.840e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.393e-06
  |ΔD|_∞ = 6.574e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5151e+00 J
  → Fracture energy : 1.6712e+00 J
  → Total energy    : 3.1863e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0116.vtu


## Step 118/200: t = 2.12e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 117 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.570e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.091e-04
  |ΔD|_∞ = 8.390e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.811e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.618e-04
  |ΔD|_∞ = 1.678e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.547e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.236e-05
  |ΔD|_∞ = 3.356e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.540e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.473e-06
  |ΔD|_∞ = 6.712e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5401e+00 J
  → Fracture energy : 1.6722e+00 J
  → Total energy    : 3.2123e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0117.vtu


## Step 119/200: t = 2.13e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 118 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.497e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.219e-04
  |ΔD|_∞ = 8.588e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.684e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.644e-04
  |ΔD|_∞ = 1.718e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.559e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.288e-05
  |ΔD|_∞ = 3.435e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.556e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.575e-06
  |ΔD|_∞ = 6.870e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5653e+00 J
  → Fracture energy : 1.6732e+00 J
  → Total energy    : 3.2384e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0118.vtu


## Step 120/200: t = 2.15e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 119 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.425e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.377e-04
  |ΔD|_∞ = 8.807e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.064e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.675e-04
  |ΔD|_∞ = 1.761e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.064e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.351e-05
  |ΔD|_∞ = 3.523e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.064e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.701e-06
  |ΔD|_∞ = 7.046e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5907e+00 J
  → Fracture energy : 1.6742e+00 J
  → Total energy    : 3.2648e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0119.vtu


## Step 121/200: t = 2.17e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 120 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.355e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.562e-04
  |ΔD|_∞ = 9.046e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.520e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.712e-04
  |ΔD|_∞ = 1.809e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.520e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.425e-05
  |ΔD|_∞ = 3.618e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.577e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.850e-06
  |ΔD|_∞ = 7.237e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6162e+00 J
  → Fracture energy : 1.6752e+00 J
  → Total energy    : 3.2914e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0120.vtu


## Step 122/200: t = 2.19e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 121 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.287e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.776e-04
  |ΔD|_∞ = 9.306e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.144e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.755e-04
  |ΔD|_∞ = 1.861e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.788e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.510e-05
  |ΔD|_∞ = 3.722e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.263e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.021e-06
  |ΔD|_∞ = 7.445e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6420e+00 J
  → Fracture energy : 1.6762e+00 J
  → Total energy    : 3.3182e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0121.vtu


## Step 123/200: t = 2.21e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 122 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.220e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.020e-04
  |ΔD|_∞ = 9.590e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.558e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.804e-04
  |ΔD|_∞ = 1.918e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.558e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.608e-05
  |ΔD|_∞ = 3.836e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.558e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.216e-06
  |ΔD|_∞ = 7.672e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6679e+00 J
  → Fracture energy : 1.6772e+00 J
  → Total energy    : 3.3451e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0122.vtu


## Step 124/200: t = 2.23e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 123 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.154e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.297e-04
  |ΔD|_∞ = 9.897e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.827e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.859e-04
  |ΔD|_∞ = 1.979e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.321e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.719e-05
  |ΔD|_∞ = 3.959e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.045e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.437e-06
  |ΔD|_∞ = 7.918e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6940e+00 J
  → Fracture energy : 1.6783e+00 J
  → Total energy    : 3.3723e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0123.vtu


## Step 125/200: t = 2.24e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 124 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.090e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.611e-04
  |ΔD|_∞ = 1.023e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.515e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.922e-04
  |ΔD|_∞ = 2.045e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.844e-05
  |ΔD|_∞ = 4.090e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.689e-06
  |ΔD|_∞ = 8.181e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7202e+00 J
  → Fracture energy : 1.6794e+00 J
  → Total energy    : 3.3997e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0124.vtu


## Step 126/200: t = 2.26e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 125 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.028e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.965e-04
  |ΔD|_∞ = 1.057e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.993e-04
  |ΔD|_∞ = 2.115e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.542e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.986e-05
  |ΔD|_∞ = 4.230e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.972e-06
  |ΔD|_∞ = 8.460e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7467e+00 J
  → Fracture energy : 1.6806e+00 J
  → Total energy    : 3.4272e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0125.vtu


## Step 127/200: t = 2.28e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 126 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.966e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.036e-03
  |ΔD|_∞ = 1.100e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.457e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.073e-04
  |ΔD|_∞ = 2.199e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.315e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.145e-05
  |ΔD|_∞ = 4.398e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.291e-06
  |ΔD|_∞ = 8.797e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7732e+00 J
  → Fracture energy : 1.6818e+00 J
  → Total energy    : 3.4550e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0126.vtu


## Step 128/200: t = 2.30e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 127 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.907e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.081e-03
  |ΔD|_∞ = 1.166e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.221e-20

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.162e-04
  |ΔD|_∞ = 2.332e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.758e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.325e-05
  |ΔD|_∞ = 4.665e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.566e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.650e-06
  |ΔD|_∞ = 9.330e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8000e+00 J
  → Fracture energy : 1.6830e+00 J
  → Total energy    : 3.4830e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0127.vtu


## Step 129/200: t = 2.32e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 128 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.848e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.132e-03
  |ΔD|_∞ = 1.239e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.085e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.263e-04
  |ΔD|_∞ = 2.479e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.526e-05
  |ΔD|_∞ = 4.957e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.296e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.052e-06
  |ΔD|_∞ = 9.914e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8269e+00 J
  → Fracture energy : 1.6843e+00 J
  → Total energy    : 3.5111e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0128.vtu


## Step 130/200: t = 2.33e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 129 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.792e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.188e-03
  |ΔD|_∞ = 1.319e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.791e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.376e-04
  |ΔD|_∞ = 2.638e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.609e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.752e-05
  |ΔD|_∞ = 5.275e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.209e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.503e-06
  |ΔD|_∞ = 1.055e-04

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8539e+00 J
  → Fracture energy : 1.6856e+00 J
  → Total energy    : 3.5395e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0129.vtu


## Step 131/200: t = 2.35e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 130 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.736e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.251e-03
  |ΔD|_∞ = 1.404e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.378e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.502e-04
  |ΔD|_∞ = 2.809e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.218e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.004e-05
  |ΔD|_∞ = 5.618e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.001e-05
  |ΔD|_∞ = 1.124e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.002e-06
  |ΔD|_∞ = 2.247e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8810e+00 J
  → Fracture energy : 1.6870e+00 J
  → Total energy    : 3.5680e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0130.vtu


## Step 132/200: t = 2.37e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 131 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.683e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.321e-03
  |ΔD|_∞ = 1.495e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.045e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.642e-04
  |ΔD|_∞ = 2.990e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.045e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.285e-05
  |ΔD|_∞ = 5.980e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.057e-05
  |ΔD|_∞ = 1.196e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.114e-06
  |ΔD|_∞ = 2.392e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9083e+00 J
  → Fracture energy : 1.6885e+00 J
  → Total energy    : 3.5968e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0131.vtu


## Step 133/200: t = 2.39e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 132 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.630e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.400e-03
  |ΔD|_∞ = 1.589e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.201e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.800e-04
  |ΔD|_∞ = 3.177e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.601e-05
  |ΔD|_∞ = 6.354e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.370e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.120e-05
  |ΔD|_∞ = 1.271e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.367e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.240e-06
  |ΔD|_∞ = 2.542e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9357e+00 J
  → Fracture energy : 1.6900e+00 J
  → Total energy    : 3.6257e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0132.vtu


## Step 134/200: t = 2.41e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 133 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.580e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.490e-03
  |ΔD|_∞ = 1.682e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.981e-04
  |ΔD|_∞ = 3.365e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.171e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.962e-05
  |ΔD|_∞ = 6.729e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.693e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.192e-05
  |ΔD|_∞ = 1.346e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.385e-06
  |ΔD|_∞ = 2.692e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9631e+00 J
  → Fracture energy : 1.6916e+00 J
  → Total energy    : 3.6548e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0133.vtu


## Step 135/200: t = 2.42e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 134 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.531e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.592e-03
  |ΔD|_∞ = 1.772e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.896e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.184e-04
  |ΔD|_∞ = 3.543e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.346e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.367e-05
  |ΔD|_∞ = 7.086e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.273e-05
  |ΔD|_∞ = 1.417e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.150e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.547e-06
  |ΔD|_∞ = 2.835e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9907e+00 J
  → Fracture energy : 1.6934e+00 J
  → Total energy    : 3.6841e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0134.vtu


## Step 136/200: t = 2.44e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 135 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.484e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.704e-03
  |ΔD|_∞ = 1.895e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.785e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.409e-04
  |ΔD|_∞ = 3.790e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.785e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.818e-05
  |ΔD|_∞ = 7.580e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.785e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.364e-05
  |ΔD|_∞ = 1.516e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.835e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.727e-06
  |ΔD|_∞ = 3.032e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0183e+00 J
  → Fracture energy : 1.6952e+00 J
  → Total energy    : 3.7135e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0135.vtu


## Step 137/200: t = 2.46e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 136 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.439e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.827e-03
  |ΔD|_∞ = 2.055e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.714e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.654e-04
  |ΔD|_∞ = 4.109e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.309e-05
  |ΔD|_∞ = 8.218e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.462e-05
  |ΔD|_∞ = 1.644e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.923e-06
  |ΔD|_∞ = 3.287e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0460e+00 J
  → Fracture energy : 1.6971e+00 J
  → Total energy    : 3.7432e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0136.vtu


## Step 138/200: t = 2.48e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 137 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.395e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.960e-03
  |ΔD|_∞ = 2.219e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.765e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.920e-04
  |ΔD|_∞ = 4.438e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.879e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.841e-05
  |ΔD|_∞ = 8.876e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.072e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.568e-05
  |ΔD|_∞ = 1.775e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.744e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.136e-06
  |ΔD|_∞ = 3.551e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0738e+00 J
  → Fracture energy : 1.6992e+00 J
  → Total energy    : 3.7730e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0137.vtu


## Step 139/200: t = 2.50e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 138 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.353e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.103e-03
  |ΔD|_∞ = 2.381e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.432e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.206e-04
  |ΔD|_∞ = 4.761e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.432e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.411e-05
  |ΔD|_∞ = 9.523e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.432e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.682e-05
  |ΔD|_∞ = 1.905e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.138e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.365e-06
  |ΔD|_∞ = 3.809e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1016e+00 J
  → Fracture energy : 1.7013e+00 J
  → Total energy    : 3.8029e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0138.vtu


## Step 140/200: t = 2.51e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 139 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.312e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.254e-03
  |ΔD|_∞ = 2.530e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.059e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.509e-04
  |ΔD|_∞ = 5.059e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.444e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.017e-05
  |ΔD|_∞ = 1.012e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.801e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.803e-05
  |ΔD|_∞ = 2.024e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.816e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.607e-06
  |ΔD|_∞ = 4.047e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1294e+00 J
  → Fracture energy : 1.7036e+00 J
  → Total energy    : 3.8331e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0139.vtu


## Step 141/200: t = 2.53e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 140 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.271e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.418e-03
  |ΔD|_∞ = 2.655e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.836e-04
  |ΔD|_∞ = 5.310e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.429e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.672e-05
  |ΔD|_∞ = 1.062e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.934e-05
  |ΔD|_∞ = 2.124e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.250e-20

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.869e-06
  |ΔD|_∞ = 4.248e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1573e+00 J
  → Fracture energy : 1.7061e+00 J
  → Total energy    : 3.8634e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0140.vtu


## Step 142/200: t = 2.55e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 141 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.232e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.603e-03
  |ΔD|_∞ = 2.743e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.714e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.205e-04
  |ΔD|_∞ = 5.485e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.714e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.041e-04
  |ΔD|_∞ = 1.097e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.082e-05
  |ΔD|_∞ = 2.194e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.724e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.164e-06
  |ΔD|_∞ = 4.388e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1852e+00 J
  → Fracture energy : 1.7086e+00 J
  → Total energy    : 3.8938e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0141.vtu


## Step 143/200: t = 2.57e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 142 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.195e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.815e-03
  |ΔD|_∞ = 2.933e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.868e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.629e-04
  |ΔD|_∞ = 5.866e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.238e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.126e-04
  |ΔD|_∞ = 1.173e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.893e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.252e-05
  |ΔD|_∞ = 2.346e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.247e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.503e-06
  |ΔD|_∞ = 4.693e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2131e+00 J
  → Fracture energy : 1.7113e+00 J
  → Total energy    : 3.9244e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0142.vtu


## Step 144/200: t = 2.59e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 143 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.162e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.063e-03
  |ΔD|_∞ = 3.152e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.314e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.126e-04
  |ΔD|_∞ = 6.303e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.132e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.225e-04
  |ΔD|_∞ = 1.261e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.872e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.450e-05
  |ΔD|_∞ = 2.521e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.314e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.901e-06
  |ΔD|_∞ = 5.043e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2409e+00 J
  → Fracture energy : 1.7142e+00 J
  → Total energy    : 3.9551e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0143.vtu


## Step 145/200: t = 2.61e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 144 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.131e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.360e-03
  |ΔD|_∞ = 3.439e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.719e-04
  |ΔD|_∞ = 6.878e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.874e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.344e-04
  |ΔD|_∞ = 1.376e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.975e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.688e-05
  |ΔD|_∞ = 2.751e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.375e-06
  |ΔD|_∞ = 5.502e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2686e+00 J
  → Fracture energy : 1.7173e+00 J
  → Total energy    : 3.9860e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0144.vtu


## Step 146/200: t = 2.62e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 145 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.106e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.710e-03
  |ΔD|_∞ = 3.795e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.716e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.419e-04
  |ΔD|_∞ = 7.590e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.716e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.484e-04
  |ΔD|_∞ = 1.518e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.716e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.968e-05
  |ΔD|_∞ = 3.036e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.935e-06
  |ΔD|_∞ = 6.072e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2962e+00 J
  → Fracture energy : 1.7206e+00 J
  → Total energy    : 4.0169e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0145.vtu


## Step 147/200: t = 2.64e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 146 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.087e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.102e-03
  |ΔD|_∞ = 4.118e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.039e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.203e-04
  |ΔD|_∞ = 8.236e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.376e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.641e-04
  |ΔD|_∞ = 1.647e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.158e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.281e-05
  |ΔD|_∞ = 3.295e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.026e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.563e-06
  |ΔD|_∞ = 6.589e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3237e+00 J
  → Fracture energy : 1.7242e+00 J
  → Total energy    : 4.0479e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0146.vtu


## Step 148/200: t = 2.66e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 147 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.069e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.546e-03
  |ΔD|_∞ = 4.296e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.375e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.091e-04
  |ΔD|_∞ = 8.591e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.395e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.818e-04
  |ΔD|_∞ = 1.718e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.395e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.637e-05
  |ΔD|_∞ = 3.437e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.413e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.273e-06
  |ΔD|_∞ = 6.873e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3510e+00 J
  → Fracture energy : 1.7280e+00 J
  → Total energy    : 4.0790e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0147.vtu


## Step 149/200: t = 2.68e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 148 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.052e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.029e-03
  |ΔD|_∞ = 4.736e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.493e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.006e-03
  |ΔD|_∞ = 9.471e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.012e-04
  |ΔD|_∞ = 1.894e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.023e-05
  |ΔD|_∞ = 3.788e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.259e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.047e-06
  |ΔD|_∞ = 7.577e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3781e+00 J
  → Fracture energy : 1.7320e+00 J
  → Total energy    : 4.1101e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0148.vtu


## Step 150/200: t = 2.70e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 149 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.035e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.529e-03
  |ΔD|_∞ = 5.186e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.072e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.106e-03
  |ΔD|_∞ = 1.037e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.072e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.212e-04
  |ΔD|_∞ = 2.074e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.063e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.424e-05
  |ΔD|_∞ = 4.149e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.065e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.847e-06
  |ΔD|_∞ = 8.297e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4051e+00 J
  → Fracture energy : 1.7362e+00 J
  → Total energy    : 4.1413e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0149.vtu


## Step 151/200: t = 2.71e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 150 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.015e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.028e-03
  |ΔD|_∞ = 5.431e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.008e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.206e-03
  |ΔD|_∞ = 1.086e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.547e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.411e-04
  |ΔD|_∞ = 2.172e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.548e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.822e-05
  |ΔD|_∞ = 4.345e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.645e-06
  |ΔD|_∞ = 8.690e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4320e+00 J
  → Fracture energy : 1.7406e+00 J
  → Total energy    : 4.1726e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0150.vtu


## Step 152/200: t = 2.73e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 151 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.994e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.531e-03
  |ΔD|_∞ = 5.457e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.241e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.306e-03
  |ΔD|_∞ = 1.091e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.612e-04
  |ΔD|_∞ = 2.183e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.412e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.225e-05
  |ΔD|_∞ = 4.366e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.419e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.045e-05
  |ΔD|_∞ = 8.732e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.479e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.090e-06
  |ΔD|_∞ = 1.746e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4587e+00 J
  → Fracture energy : 1.7452e+00 J
  → Total energy    : 4.2039e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0151.vtu


## Step 153/200: t = 2.75e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 152 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.975e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.070e-03
  |ΔD|_∞ = 5.766e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.093e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.414e-03
  |ΔD|_∞ = 1.153e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.828e-04
  |ΔD|_∞ = 2.306e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.353e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.656e-05
  |ΔD|_∞ = 4.613e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.353e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.131e-05
  |ΔD|_∞ = 9.226e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.263e-06
  |ΔD|_∞ = 1.845e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4853e+00 J
  → Fracture energy : 1.7499e+00 J
  → Total energy    : 4.2352e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0152.vtu


## Step 154/200: t = 2.77e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 153 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.959e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.593e-03
  |ΔD|_∞ = 5.784e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.540e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.519e-03
  |ΔD|_∞ = 1.157e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.037e-04
  |ΔD|_∞ = 2.314e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.424e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.074e-05
  |ΔD|_∞ = 4.628e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.424e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.215e-05
  |ΔD|_∞ = 9.255e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.166e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.430e-06
  |ΔD|_∞ = 1.851e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5119e+00 J
  → Fracture energy : 1.7548e+00 J
  → Total energy    : 4.2666e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0153.vtu


## Step 155/200: t = 2.79e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 154 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.937e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.138e-03
  |ΔD|_∞ = 5.985e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.628e-03
  |ΔD|_∞ = 1.197e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.255e-04
  |ΔD|_∞ = 2.394e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.511e-05
  |ΔD|_∞ = 4.788e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.373e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.302e-05
  |ΔD|_∞ = 9.575e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.424e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.604e-06
  |ΔD|_∞ = 1.915e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5383e+00 J
  → Fracture energy : 1.7598e+00 J
  → Total energy    : 4.2981e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0154.vtu


## Step 156/200: t = 2.80e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 155 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.920e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.775e-03
  |ΔD|_∞ = 6.084e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.755e-03
  |ΔD|_∞ = 1.217e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.119e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.510e-04
  |ΔD|_∞ = 2.434e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.020e-05
  |ΔD|_∞ = 4.867e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.119e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.404e-05
  |ΔD|_∞ = 9.735e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.119e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.808e-06
  |ΔD|_∞ = 1.947e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5646e+00 J
  → Fracture energy : 1.7650e+00 J
  → Total energy    : 4.3295e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0155.vtu


## Step 157/200: t = 2.82e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 156 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.912e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.444e-03
  |ΔD|_∞ = 6.319e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.116e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.889e-03
  |ΔD|_∞ = 1.264e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.117e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.778e-04
  |ΔD|_∞ = 2.528e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.125e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.555e-05
  |ΔD|_∞ = 5.056e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.270e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.511e-05
  |ΔD|_∞ = 1.011e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.022e-06
  |ΔD|_∞ = 2.022e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5907e+00 J
  → Fracture energy : 1.7703e+00 J
  → Total energy    : 4.3610e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0156.vtu


## Step 158/200: t = 2.84e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 157 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.890e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.010e-02
  |ΔD|_∞ = 6.415e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.222e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.021e-03
  |ΔD|_∞ = 1.283e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.623e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.042e-04
  |ΔD|_∞ = 2.566e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.215e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.084e-05
  |ΔD|_∞ = 5.132e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.961e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.617e-05
  |ΔD|_∞ = 1.026e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.080e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.233e-06
  |ΔD|_∞ = 2.053e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6168e+00 J
  → Fracture energy : 1.7757e+00 J
  → Total energy    : 4.3925e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0157.vtu


## Step 159/200: t = 2.86e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 158 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.892e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.087e-02
  |ΔD|_∞ = 6.732e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.542e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.175e-03
  |ΔD|_∞ = 1.346e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.542e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.350e-04
  |ΔD|_∞ = 2.693e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.700e-05
  |ΔD|_∞ = 5.385e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.542e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.740e-05
  |ΔD|_∞ = 1.077e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.463e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.480e-06
  |ΔD|_∞ = 2.154e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6426e+00 J
  → Fracture energy : 1.7814e+00 J
  → Total energy    : 4.4240e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0158.vtu


## Step 160/200: t = 2.88e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 159 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.899e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.169e-02
  |ΔD|_∞ = 6.926e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.615e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.338e-03
  |ΔD|_∞ = 1.385e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.675e-04
  |ΔD|_∞ = 2.770e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.127e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.351e-05
  |ΔD|_∞ = 5.541e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.870e-05
  |ΔD|_∞ = 1.108e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.740e-06
  |ΔD|_∞ = 2.216e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6682e+00 J
  → Fracture energy : 1.7872e+00 J
  → Total energy    : 4.4554e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0159.vtu


## Step 161/200: t = 2.89e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 160 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.887e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.237e-02
  |ΔD|_∞ = 7.177e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.629e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.475e-03
  |ΔD|_∞ = 1.435e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.629e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.950e-04
  |ΔD|_∞ = 2.871e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.144e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.900e-05
  |ΔD|_∞ = 5.741e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.602e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.980e-05
  |ΔD|_∞ = 1.148e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.623e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.960e-06
  |ΔD|_∞ = 2.297e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6937e+00 J
  → Fracture energy : 1.7932e+00 J
  → Total energy    : 4.4869e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0160.vtu


## Step 162/200: t = 2.91e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 161 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.892e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.294e-02
  |ΔD|_∞ = 7.338e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.286e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.587e-03
  |ΔD|_∞ = 1.468e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.668e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.174e-04
  |ΔD|_∞ = 2.935e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.404e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.035e-04
  |ΔD|_∞ = 5.870e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.404e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.070e-05
  |ΔD|_∞ = 1.174e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.668e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.139e-06
  |ΔD|_∞ = 2.348e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7190e+00 J
  → Fracture energy : 1.7993e+00 J
  → Total energy    : 4.5183e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0161.vtu


## Step 163/200: t = 2.93e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 162 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.908e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.323e-02
  |ΔD|_∞ = 7.453e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.541e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.646e-03
  |ΔD|_∞ = 1.491e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.541e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.292e-04
  |ΔD|_∞ = 2.981e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.541e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.058e-04
  |ΔD|_∞ = 5.963e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.658e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.117e-05
  |ΔD|_∞ = 1.193e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.658e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.234e-06
  |ΔD|_∞ = 2.385e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7441e+00 J
  → Fracture energy : 1.8056e+00 J
  → Total energy    : 4.5496e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0162.vtu


## Step 164/200: t = 2.95e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 163 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.923e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.323e-02
  |ΔD|_∞ = 7.584e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.215e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.646e-03
  |ΔD|_∞ = 1.517e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.121e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.292e-04
  |ΔD|_∞ = 3.033e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.601e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.058e-04
  |ΔD|_∞ = 6.067e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.707e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.117e-05
  |ΔD|_∞ = 1.213e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.969e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.233e-06
  |ΔD|_∞ = 2.427e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7689e+00 J
  → Fracture energy : 1.8120e+00 J
  → Total energy    : 4.5810e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0163.vtu


## Step 165/200: t = 2.97e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 164 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.943e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.293e-02
  |ΔD|_∞ = 7.736e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.678e-20

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.586e-03
  |ΔD|_∞ = 1.547e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.171e-04
  |ΔD|_∞ = 3.094e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.482e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.034e-04
  |ΔD|_∞ = 6.189e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.482e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.068e-05
  |ΔD|_∞ = 1.238e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.137e-06
  |ΔD|_∞ = 2.475e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7936e+00 J
  → Fracture energy : 1.8187e+00 J
  → Total energy    : 4.6123e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0164.vtu


## Step 166/200: t = 2.98e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 165 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.991e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.233e-02
  |ΔD|_∞ = 7.926e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.086e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.467e-03
  |ΔD|_∞ = 1.585e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.086e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.933e-04
  |ΔD|_∞ = 3.170e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.083e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.866e-05
  |ΔD|_∞ = 6.340e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.122e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.973e-05
  |ΔD|_∞ = 1.268e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.919e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.947e-06
  |ΔD|_∞ = 2.536e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8182e+00 J
  → Fracture energy : 1.8255e+00 J
  → Total energy    : 4.6436e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0165.vtu


## Step 167/200: t = 3.00e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 166 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.078e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.168e-02
  |ΔD|_∞ = 8.011e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.344e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.336e-03
  |ΔD|_∞ = 1.602e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.077e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.672e-04
  |ΔD|_∞ = 3.205e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.345e-05
  |ΔD|_∞ = 6.409e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.869e-05
  |ΔD|_∞ = 1.282e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.762e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.738e-06
  |ΔD|_∞ = 2.564e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8425e+00 J
  → Fracture energy : 1.8324e+00 J
  → Total energy    : 4.6749e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0166.vtu


## Step 168/200: t = 3.02e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 167 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.116e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.116e-02
  |ΔD|_∞ = 8.089e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.043e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.232e-03
  |ΔD|_∞ = 1.618e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.465e-04
  |ΔD|_∞ = 3.236e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.930e-05
  |ΔD|_∞ = 6.471e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.504e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.786e-05
  |ΔD|_∞ = 1.294e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.504e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.572e-06
  |ΔD|_∞ = 2.589e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8666e+00 J
  → Fracture energy : 1.8395e+00 J
  → Total energy    : 4.7061e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0167.vtu


## Step 169/200: t = 3.04e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 168 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.038e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.078e-02
  |ΔD|_∞ = 8.117e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.817e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.155e-03
  |ΔD|_∞ = 1.623e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.000e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.310e-04
  |ΔD|_∞ = 3.247e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.031e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.619e-05
  |ΔD|_∞ = 6.494e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.047e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.724e-05
  |ΔD|_∞ = 1.299e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.396e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.448e-06
  |ΔD|_∞ = 2.598e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8905e+00 J
  → Fracture energy : 1.8469e+00 J
  → Total energy    : 4.7374e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0168.vtu


## Step 170/200: t = 3.06e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 169 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.957e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.054e-02
  |ΔD|_∞ = 8.116e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.594e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.108e-03
  |ΔD|_∞ = 1.623e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.217e-04
  |ΔD|_∞ = 3.246e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.925e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.433e-05
  |ΔD|_∞ = 6.493e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.888e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.687e-05
  |ΔD|_∞ = 1.299e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.773e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.373e-06
  |ΔD|_∞ = 2.597e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9142e+00 J
  → Fracture energy : 1.8544e+00 J
  → Total energy    : 4.7686e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0169.vtu


## Step 171/200: t = 3.08e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 170 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.938e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.058e-02
  |ΔD|_∞ = 8.057e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.116e-03
  |ΔD|_∞ = 1.611e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.957e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.231e-04
  |ΔD|_∞ = 3.223e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.957e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.463e-05
  |ΔD|_∞ = 6.446e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.693e-05
  |ΔD|_∞ = 1.289e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.186e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.385e-06
  |ΔD|_∞ = 2.578e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9377e+00 J
  → Fracture energy : 1.8621e+00 J
  → Total energy    : 4.7998e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0170.vtu


## Step 172/200: t = 3.09e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 171 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.894e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.091e-02
  |ΔD|_∞ = 8.104e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.078e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.181e-03
  |ΔD|_∞ = 1.621e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.189e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.361e-04
  |ΔD|_∞ = 3.242e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.112e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.722e-05
  |ΔD|_∞ = 6.483e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.189e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.744e-05
  |ΔD|_∞ = 1.297e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.489e-06
  |ΔD|_∞ = 2.593e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9609e+00 J
  → Fracture energy : 1.8700e+00 J
  → Total energy    : 4.8309e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0171.vtu


## Step 173/200: t = 3.11e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 172 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.839e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.163e-02
  |ΔD|_∞ = 8.202e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.193e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.326e-03
  |ΔD|_∞ = 1.640e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.542e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.651e-04
  |ΔD|_∞ = 3.281e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.970e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.303e-05
  |ΔD|_∞ = 6.561e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.970e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.861e-05
  |ΔD|_∞ = 1.312e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.519e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.721e-06
  |ΔD|_∞ = 2.625e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9840e+00 J
  → Fracture energy : 1.8781e+00 J
  → Total energy    : 4.8621e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0172.vtu


## Step 174/200: t = 3.13e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 173 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.833e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.260e-02
  |ΔD|_∞ = 8.481e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.519e-03
  |ΔD|_∞ = 1.696e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.604e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.038e-04
  |ΔD|_∞ = 3.392e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.600e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.008e-04
  |ΔD|_∞ = 6.785e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.784e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.015e-05
  |ΔD|_∞ = 1.357e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.312e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.031e-06
  |ΔD|_∞ = 2.714e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.0068e+00 J
  → Fracture energy : 1.8864e+00 J
  → Total energy    : 4.8932e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0173.vtu


## Step 175/200: t = 3.15e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 174 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.817e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.381e-02
  |ΔD|_∞ = 8.552e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.077e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.763e-03
  |ΔD|_∞ = 1.710e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.077e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.526e-04
  |ΔD|_∞ = 3.421e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.629e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.105e-04
  |ΔD|_∞ = 6.842e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.982e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.210e-05
  |ΔD|_∞ = 1.368e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.415e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.420e-06
  |ΔD|_∞ = 2.737e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.0290e+00 J
  → Fracture energy : 1.8950e+00 J
  → Total energy    : 4.9240e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0174.vtu


## Step 176/200: t = 3.17e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 175 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.814e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.510e-02
  |ΔD|_∞ = 8.880e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.038e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.021e-03
  |ΔD|_∞ = 1.776e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.228e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.041e-04
  |ΔD|_∞ = 3.552e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.192e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.208e-04
  |ΔD|_∞ = 7.104e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.121e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.416e-05
  |ΔD|_∞ = 1.421e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.04e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.229e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.833e-06
  |ΔD|_∞ = 2.841e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.0509e+00 J
  → Fracture energy : 1.9038e+00 J
  → Total energy    : 4.9546e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0175.vtu


## Step 177/200: t = 3.18e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 176 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.801e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.620e-02
  |ΔD|_∞ = 9.075e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.502e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.240e-03
  |ΔD|_∞ = 1.815e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.345e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.479e-04
  |ΔD|_∞ = 3.630e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.345e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.296e-04
  |ΔD|_∞ = 7.260e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.345e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.592e-05
  |ΔD|_∞ = 1.452e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.08e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.345e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.184e-06
  |ΔD|_∞ = 2.904e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.0725e+00 J
  → Fracture energy : 1.9128e+00 J
  → Total energy    : 4.9852e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0176.vtu


## Step 178/200: t = 3.20e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 177 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.840e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.701e-02
  |ΔD|_∞ = 9.199e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.953e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.402e-03
  |ΔD|_∞ = 1.840e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.804e-04
  |ΔD|_∞ = 3.680e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.361e-04
  |ΔD|_∞ = 7.360e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.997e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.722e-05
  |ΔD|_∞ = 1.472e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.12e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.000e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.443e-06
  |ΔD|_∞ = 2.944e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.0938e+00 J
  → Fracture energy : 1.9219e+00 J
  → Total energy    : 5.0156e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0177.vtu


## Step 179/200: t = 3.22e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 178 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.901e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.741e-02
  |ΔD|_∞ = 9.353e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.549e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.482e-03
  |ΔD|_∞ = 1.871e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.964e-04
  |ΔD|_∞ = 3.741e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.742e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.393e-04
  |ΔD|_∞ = 7.482e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.786e-05
  |ΔD|_∞ = 1.496e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.16e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.742e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.571e-06
  |ΔD|_∞ = 2.993e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1147e+00 J
  → Fracture energy : 1.9312e+00 J
  → Total energy    : 5.0458e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0178.vtu


## Step 180/200: t = 3.24e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 179 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.958e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.722e-02
  |ΔD|_∞ = 9.465e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.823e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.445e-03
  |ΔD|_∞ = 1.893e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.823e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.888e-04
  |ΔD|_∞ = 3.786e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.823e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.378e-04
  |ΔD|_∞ = 7.572e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.936e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.755e-05
  |ΔD|_∞ = 1.514e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.939e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.510e-06
  |ΔD|_∞ = 3.029e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1352e+00 J
  → Fracture energy : 1.9407e+00 J
  → Total energy    : 5.0759e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0179.vtu


## Step 181/200: t = 3.26e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 180 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.047e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.640e-02
  |ΔD|_∞ = 9.617e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.307e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.281e-03
  |ΔD|_∞ = 1.923e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.542e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.561e-04
  |ΔD|_∞ = 3.847e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.542e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.312e-04
  |ΔD|_∞ = 7.693e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.624e-05
  |ΔD|_∞ = 1.539e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.24e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.942e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.249e-06
  |ΔD|_∞ = 3.077e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1554e+00 J
  → Fracture energy : 1.9504e+00 J
  → Total energy    : 5.1058e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0180.vtu


## Step 182/200: t = 3.27e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 181 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.195e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.523e-02
  |ΔD|_∞ = 9.669e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.291e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.047e-03
  |ΔD|_∞ = 1.934e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.498e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.094e-04
  |ΔD|_∞ = 3.868e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.857e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.219e-04
  |ΔD|_∞ = 7.735e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.857e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.437e-05
  |ΔD|_∞ = 1.547e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.28e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.952e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.875e-06
  |ΔD|_∞ = 3.094e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1753e+00 J
  → Fracture energy : 1.9603e+00 J
  → Total energy    : 5.1356e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0181.vtu


## Step 183/200: t = 3.29e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 182 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.366e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.424e-02
  |ΔD|_∞ = 9.770e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.435e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.847e-03
  |ΔD|_∞ = 1.954e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.695e-04
  |ΔD|_∞ = 3.908e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.139e-04
  |ΔD|_∞ = 7.816e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.278e-05
  |ΔD|_∞ = 1.563e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.32e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.426e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.556e-06
  |ΔD|_∞ = 3.126e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1947e+00 J
  → Fracture energy : 1.9704e+00 J
  → Total energy    : 5.1652e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0182.vtu


## Step 184/200: t = 3.31e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 183 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.319e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.369e-02
  |ΔD|_∞ = 9.830e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.831e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.738e-03
  |ΔD|_∞ = 1.966e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.426e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.475e-04
  |ΔD|_∞ = 3.932e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.298e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.095e-04
  |ΔD|_∞ = 7.864e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.184e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.190e-05
  |ΔD|_∞ = 1.573e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.36e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.300e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.380e-06
  |ΔD|_∞ = 3.146e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.2139e+00 J
  → Fracture energy : 1.9808e+00 J
  → Total energy    : 5.1947e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0183.vtu


## Step 185/200: t = 3.33e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 184 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.142e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.365e-02
  |ΔD|_∞ = 9.858e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.212e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.730e-03
  |ΔD|_∞ = 1.972e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.212e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.460e-04
  |ΔD|_∞ = 3.943e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.092e-04
  |ΔD|_∞ = 7.887e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.184e-05
  |ΔD|_∞ = 1.577e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.367e-06
  |ΔD|_∞ = 3.155e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.2327e+00 J
  → Fracture energy : 1.9914e+00 J
  → Total energy    : 5.2241e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0184.vtu


## Step 186/200: t = 3.35e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 185 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.087e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.416e-02
  |ΔD|_∞ = 1.003e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.832e-03
  |ΔD|_∞ = 2.006e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.612e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.664e-04
  |ΔD|_∞ = 4.012e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.357e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.133e-04
  |ΔD|_∞ = 8.024e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.265e-05
  |ΔD|_∞ = 1.605e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.44e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.357e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.531e-06
  |ΔD|_∞ = 3.210e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.2510e+00 J
  → Fracture energy : 2.0023e+00 J
  → Total energy    : 5.2533e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0185.vtu


## Step 187/200: t = 3.36e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 186 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.067e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.514e-02
  |ΔD|_∞ = 1.003e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.074e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.029e-03
  |ΔD|_∞ = 2.005e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.057e-04
  |ΔD|_∞ = 4.010e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.211e-04
  |ΔD|_∞ = 8.021e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.423e-05
  |ΔD|_∞ = 1.604e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.48e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.846e-06
  |ΔD|_∞ = 3.208e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.2689e+00 J
  → Fracture energy : 2.0134e+00 J
  → Total energy    : 5.2823e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0186.vtu


## Step 188/200: t = 3.38e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 187 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.015e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.645e-02
  |ΔD|_∞ = 1.009e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.139e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.290e-03
  |ΔD|_∞ = 2.018e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.581e-04
  |ΔD|_∞ = 4.036e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.217e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.316e-04
  |ΔD|_∞ = 8.072e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.217e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.632e-05
  |ΔD|_∞ = 1.614e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.52e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.063e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.265e-06
  |ΔD|_∞ = 3.229e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.2863e+00 J
  → Fracture energy : 2.0247e+00 J
  → Total energy    : 5.3110e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0187.vtu


## Step 189/200: t = 3.40e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 188 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.956e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.790e-02
  |ΔD|_∞ = 1.030e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.579e-03
  |ΔD|_∞ = 2.059e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.158e-04
  |ΔD|_∞ = 4.118e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.967e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.432e-04
  |ΔD|_∞ = 8.236e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.248e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.863e-05
  |ΔD|_∞ = 1.647e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.56e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.593e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.727e-06
  |ΔD|_∞ = 3.295e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3032e+00 J
  → Fracture energy : 2.0362e+00 J
  → Total energy    : 5.3395e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0188.vtu


## Step 190/200: t = 3.42e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 189 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.963e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.918e-02
  |ΔD|_∞ = 1.032e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.233e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.837e-03
  |ΔD|_∞ = 2.064e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.233e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.674e-04
  |ΔD|_∞ = 4.128e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.981e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.535e-04
  |ΔD|_∞ = 8.256e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.694e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.069e-05
  |ΔD|_∞ = 1.651e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.982e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.139e-06
  |ΔD|_∞ = 3.302e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3197e+00 J
  → Fracture energy : 2.0480e+00 J
  → Total energy    : 5.3677e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0189.vtu


## Step 191/200: t = 3.44e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 190 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.998e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.001e-02
  |ΔD|_∞ = 1.053e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.113e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.001e-03
  |ΔD|_∞ = 2.105e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.028e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.003e-04
  |ΔD|_∞ = 4.211e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.480e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.601e-04
  |ΔD|_∞ = 8.422e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.301e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.201e-05
  |ΔD|_∞ = 1.684e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.64e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.912e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.402e-06
  |ΔD|_∞ = 3.369e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3358e+00 J
  → Fracture energy : 2.0600e+00 J
  → Total energy    : 5.3957e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0190.vtu


## Step 192/200: t = 3.46e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 191 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.070e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.003e-02
  |ΔD|_∞ = 1.059e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.983e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.007e-03
  |ΔD|_∞ = 2.117e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.983e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.014e-04
  |ΔD|_∞ = 4.234e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.981e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.603e-04
  |ΔD|_∞ = 8.468e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.206e-05
  |ΔD|_∞ = 1.694e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.68e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.086e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.411e-06
  |ΔD|_∞ = 3.387e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3514e+00 J
  → Fracture energy : 2.0721e+00 J
  → Total energy    : 5.4235e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0191.vtu


## Step 193/200: t = 3.47e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 192 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.192e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.909e-02
  |ΔD|_∞ = 1.066e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.696e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.818e-03
  |ΔD|_∞ = 2.131e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.696e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.637e-04
  |ΔD|_∞ = 4.263e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.527e-04
  |ΔD|_∞ = 8.526e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.271e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.055e-05
  |ΔD|_∞ = 1.705e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.72e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.270e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.109e-06
  |ΔD|_∞ = 3.410e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3665e+00 J
  → Fracture energy : 2.0845e+00 J
  → Total energy    : 5.4510e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0192.vtu


## Step 194/200: t = 3.49e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 193 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.423e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.747e-02
  |ΔD|_∞ = 1.078e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.493e-03
  |ΔD|_∞ = 2.155e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.251e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.987e-04
  |ΔD|_∞ = 4.311e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.251e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.397e-04
  |ΔD|_∞ = 8.622e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.086e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.795e-05
  |ΔD|_∞ = 1.724e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.76e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.589e-06
  |ΔD|_∞ = 3.449e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3812e+00 J
  → Fracture energy : 2.0971e+00 J
  → Total energy    : 5.4783e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0193.vtu


## Step 195/200: t = 3.51e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 194 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.517e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.605e-02
  |ΔD|_∞ = 1.082e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.526e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.211e-03
  |ΔD|_∞ = 2.164e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.881e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.422e-04
  |ΔD|_∞ = 4.327e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.881e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.284e-04
  |ΔD|_∞ = 8.654e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.881e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.569e-05
  |ΔD|_∞ = 1.731e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.138e-06
  |ΔD|_∞ = 3.462e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3954e+00 J
  → Fracture energy : 2.1099e+00 J
  → Total energy    : 5.5053e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0194.vtu


## Step 196/200: t = 3.53e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 195 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.651e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.558e-02
  |ΔD|_∞ = 1.081e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.057e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.115e-03
  |ΔD|_∞ = 2.162e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.230e-04
  |ΔD|_∞ = 4.324e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.734e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.246e-04
  |ΔD|_∞ = 8.648e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.684e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.492e-05
  |ΔD|_∞ = 1.730e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.84e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.146e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.984e-06
  |ΔD|_∞ = 3.459e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.4091e+00 J
  → Fracture energy : 2.1231e+00 J
  → Total energy    : 5.5322e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0195.vtu


## Step 197/200: t = 3.55e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 196 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.445e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.610e-02
  |ΔD|_∞ = 1.078e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.360e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.220e-03
  |ΔD|_∞ = 2.156e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.098e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.439e-04
  |ΔD|_∞ = 4.312e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.288e-04
  |ΔD|_∞ = 8.624e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.098e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.576e-05
  |ΔD|_∞ = 1.725e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.88e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.090e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.151e-06
  |ΔD|_∞ = 3.450e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.4224e+00 J
  → Fracture energy : 2.1365e+00 J
  → Total energy    : 5.5589e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0196.vtu


## Step 198/200: t = 3.56e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 197 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.334e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.738e-02
  |ΔD|_∞ = 1.098e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.980e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.476e-03
  |ΔD|_∞ = 2.196e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.844e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.951e-04
  |ΔD|_∞ = 4.391e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.894e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.390e-04
  |ΔD|_∞ = 8.783e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.110e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.780e-05
  |ΔD|_∞ = 1.757e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.92e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.894e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.561e-06
  |ΔD|_∞ = 3.513e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.4352e+00 J
  → Fracture energy : 2.1501e+00 J
  → Total energy    : 5.5853e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0197.vtu


## Step 199/200: t = 3.58e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 198 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.213e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.904e-02
  |ΔD|_∞ = 1.111e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.747e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.808e-03
  |ΔD|_∞ = 2.223e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.748e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.615e-04
  |ΔD|_∞ = 4.445e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.523e-04
  |ΔD|_∞ = 8.891e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.831e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.046e-05
  |ΔD|_∞ = 1.778e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.96e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.831e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.092e-06
  |ΔD|_∞ = 3.556e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.4472e+00 J
  → Fracture energy : 2.1641e+00 J
  → Total energy    : 5.6113e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0198.vtu


## Step 200/200: t = 3.60e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 199 | dt = 1.81e+01 s
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
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.131e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.063e-02
  |ΔD|_∞ = 1.131e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.020e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.126e-03
  |ΔD|_∞ = 2.262e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.193e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.252e-04
  |ΔD|_∞ = 4.525e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.019e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.650e-04
  |ΔD|_∞ = 9.050e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.301e-05
  |ΔD|_∞ = 1.810e-04

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.019e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.602e-06
  |ΔD|_∞ = 3.620e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.4586e+00 J
  → Fracture energy : 2.1783e+00 J
  → Total energy    : 5.6369e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0199.vtu

Simulation completed in 571.99 s
Total time steps solved: 200
