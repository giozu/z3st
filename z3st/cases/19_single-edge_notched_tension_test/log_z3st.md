Info    : Reading 'mesh.msh'...
Info    : 60634 nodes
Info    : 121891 elements
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
  → Time steps          : 141
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
  **[INFO]** Step-dependent Dirichlet list (2D), length 141
  **[INFO]** Dirichlet mechanical BC on 'steel' → [0.0, 0.0] at region 'ymax'

Setting damage boundary conditions...
  **[INFO]** Dirichlet damage BC on 'steel' → D = 1.0 at region 'crack'


## Step 01/141: t = 0.00e+00 s | LHR = 0.00e+00 W/m

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
  → Fracture energy : 1.3581e+00 J
  → Total energy    : 1.3581e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0000.vtu


## Step 02/141: t = 2.57e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 1 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.169e-05
  |ΔD|_∞ = 1.326e-05

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5e-08]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.324e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.339e-06
  |ΔD|_∞ = 2.652e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7981e-04 J
  → Fracture energy : 1.3581e+00 J
  → Total energy    : 1.3584e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0001.vtu


## Step 03/141: t = 5.14e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 2 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.000e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.448e-05
  |ΔD|_∞ = 3.929e-05

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.726e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.897e-06
  |ΔD|_∞ = 7.859e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1192e-03 J
  → Fracture energy : 1.3581e+00 J
  → Total energy    : 1.3592e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0002.vtu


## Step 04/141: t = 7.71e+01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 3 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.333e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.810e-05
  |ΔD|_∞ = 6.615e-05

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.5e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.715e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.162e-05
  |ΔD|_∞ = 1.323e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.5e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.597e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.324e-06
  |ΔD|_∞ = 2.646e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5181e-03 J
  → Fracture energy : 1.3581e+00 J
  → Total energy    : 1.3606e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0003.vtu


## Step 05/141: t = 1.03e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 4 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  ||Δu||/||u|| = 2.500e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.995e-05
  |ΔD|_∞ = 9.117e-05

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.595e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.599e-05
  |ΔD|_∞ = 1.823e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.198e-06
  |ΔD|_∞ = 3.647e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.4763e-03 J
  → Fracture energy : 1.3581e+00 J
  → Total energy    : 1.3626e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0004.vtu


## Step 06/141: t = 1.29e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 5 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.000e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.030e-04
  |ΔD|_∞ = 1.175e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.5e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.060e-05
  |ΔD|_∞ = 2.350e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.5e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.639e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.120e-06
  |ΔD|_∞ = 4.700e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.9936e-03 J
  → Fracture energy : 1.3581e+00 J
  → Total energy    : 1.3651e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0005.vtu


## Step 07/141: t = 1.54e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 6 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.667e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.262e-04
  |ΔD|_∞ = 1.441e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.573e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.524e-05
  |ΔD|_∞ = 2.882e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.574e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.047e-06
  |ΔD|_∞ = 5.763e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0070e-02 J
  → Fracture energy : 1.3581e+00 J
  → Total energy    : 1.3682e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0006.vtu


## Step 08/141: t = 1.80e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 7 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.429e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.495e-04
  |ΔD|_∞ = 1.709e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.5e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.061e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.991e-05
  |ΔD|_∞ = 3.418e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.5e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.701e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.982e-06
  |ΔD|_∞ = 6.837e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3704e-02 J
  → Fracture energy : 1.3581e+00 J
  → Total energy    : 1.3718e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0007.vtu


## Step 09/141: t = 2.06e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 8 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  ||Δu||/||u|| = 1.250e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.731e-04
  |ΔD|_∞ = 1.981e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.466e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.462e-05
  |ΔD|_∞ = 3.962e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.925e-06
  |ΔD|_∞ = 7.924e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7896e-02 J
  → Fracture energy : 1.3581e+00 J
  → Total energy    : 1.3760e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0008.vtu


## Step 10/141: t = 2.31e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 9 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.111e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.970e-04
  |ΔD|_∞ = 2.257e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.5e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.788e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.940e-05
  |ΔD|_∞ = 4.514e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.5e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.788e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.879e-06
  |ΔD|_∞ = 9.028e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2646e-02 J
  → Fracture energy : 1.3581e+00 J
  → Total energy    : 1.3807e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0009.vtu


## Step 11/141: t = 2.57e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 10 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.000e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.211e-04
  |ΔD|_∞ = 2.538e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.684e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.422e-05
  |ΔD|_∞ = 5.075e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.845e-06
  |ΔD|_∞ = 1.015e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7952e-02 J
  → Fracture energy : 1.3581e+00 J
  → Total energy    : 1.3861e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0010.vtu


## Step 12/141: t = 2.83e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 11 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.095e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.455e-04
  |ΔD|_∞ = 2.826e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.5e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.910e-05
  |ΔD|_∞ = 5.651e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.5e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.820e-06
  |ΔD|_∞ = 1.130e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 3 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3815e-02 J
  → Fracture energy : 1.3581e+00 J
  → Total energy    : 1.3919e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0011.vtu


## Step 13/141: t = 3.09e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 12 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  ||Δu||/||u|| = 8.338e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.705e-04
  |ΔD|_∞ = 3.118e-04

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
  ||ΔD||/||D|| = 5.411e-05
  |ΔD|_∞ = 6.235e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.082e-05
  |ΔD|_∞ = 1.247e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.987e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.164e-06
  |ΔD|_∞ = 2.494e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.0232e-02 J
  → Fracture energy : 1.3581e+00 J
  → Total energy    : 1.3983e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0012.vtu


## Step 14/141: t = 3.34e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 13 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.697e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.942e-04
  |ΔD|_∞ = 3.397e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.5e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.884e-05
  |ΔD|_∞ = 6.794e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.5e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.177e-05
  |ΔD|_∞ = 1.359e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 6.5e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.353e-06
  |ΔD|_∞ = 2.718e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.7204e-02 J
  → Fracture energy : 1.3581e+00 J
  → Total energy    : 1.4053e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0013.vtu


## Step 15/141: t = 3.60e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 14 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.148e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.199e-04
  |ΔD|_∞ = 3.702e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.847e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.398e-05
  |ΔD|_∞ = 7.404e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.376e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.280e-05
  |ΔD|_∞ = 1.481e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.458e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.559e-06
  |ΔD|_∞ = 2.962e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.4729e-02 J
  → Fracture energy : 1.3581e+00 J
  → Total energy    : 1.4129e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0014.vtu


## Step 16/141: t = 3.86e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 15 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.672e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.462e-04
  |ΔD|_∞ = 4.016e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.5e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.715e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.924e-05
  |ΔD|_∞ = 8.033e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.5e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.224e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.385e-05
  |ΔD|_∞ = 1.607e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 7.5e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.616e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.769e-06
  |ΔD|_∞ = 3.213e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.2806e-02 J
  → Fracture energy : 1.3582e+00 J
  → Total energy    : 1.4210e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0015.vtu


## Step 17/141: t = 4.11e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 16 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  ||Δu||/||u|| = 6.256e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.731e-04
  |ΔD|_∞ = 4.342e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 8e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.462e-05
  |ΔD|_∞ = 8.683e-05

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
  ||ΔD||/||D|| = 1.492e-05
  |ΔD|_∞ = 1.737e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 8e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.985e-06
  |ΔD|_∞ = 3.473e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.1434e-02 J
  → Fracture energy : 1.3582e+00 J
  → Total energy    : 1.4296e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0016.vtu


## Step 18/141: t = 4.37e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 17 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.889e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.008e-04
  |ΔD|_∞ = 4.678e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 8.5e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.798e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.015e-05
  |ΔD|_∞ = 9.356e-05

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 8.5e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.603e-05
  |ΔD|_∞ = 1.871e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 8.5e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
