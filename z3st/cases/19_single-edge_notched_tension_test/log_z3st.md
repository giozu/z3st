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
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.206e-06
  |ΔD|_∞ = 3.742e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.0611e-02 J
  → Fracture energy : 1.3582e+00 J
  → Total energy    : 1.4388e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0017.vtu


## Step 19/141: t = 4.63e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 18 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.563e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.292e-04
  |ΔD|_∞ = 5.027e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 9e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.460e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.584e-05
  |ΔD|_∞ = 1.005e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 9e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.460e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.717e-05
  |ΔD|_∞ = 2.011e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 9e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.434e-06
  |ΔD|_∞ = 4.021e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.0336e-02 J
  → Fracture energy : 1.3582e+00 J
  → Total energy    : 1.4486e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0018.vtu


## Step 20/141: t = 4.89e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 19 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.272e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.586e-04
  |ΔD|_∞ = 5.390e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 9.5e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.172e-05
  |ΔD|_∞ = 1.078e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 9.5e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.752e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.834e-05
  |ΔD|_∞ = 2.156e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 9.5e-07]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.669e-06
  |ΔD|_∞ = 4.312e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0061e-01 J
  → Fracture energy : 1.3583e+00 J
  → Total energy    : 1.4589e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0019.vtu


## Step 21/141: t = 5.14e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 20 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  ||Δu||/||u|| = 5.009e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.890e-04
  |ΔD|_∞ = 5.793e-04

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
  ||ΔD||/||D|| = 9.779e-05
  |ΔD|_∞ = 1.159e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.956e-05
  |ΔD|_∞ = 2.317e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.912e-06
  |ΔD|_∞ = 4.634e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1142e-01 J
  → Fracture energy : 1.3583e+00 J
  → Total energy    : 1.4698e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0020.vtu


## Step 22/141: t = 5.40e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 21 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.05e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.772e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.206e-04
  |ΔD|_∞ = 6.224e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.05e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.132e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.041e-04
  |ΔD|_∞ = 1.245e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.05e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.132e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.082e-05
  |ΔD|_∞ = 2.489e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.05e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.621e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.165e-06
  |ΔD|_∞ = 4.979e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2278e-01 J
  → Fracture energy : 1.3584e+00 J
  → Total energy    : 1.4812e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0021.vtu


## Step 23/141: t = 5.66e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 22 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.1e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.556e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.533e-04
  |ΔD|_∞ = 6.674e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.1e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.491e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.107e-04
  |ΔD|_∞ = 1.335e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.1e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.708e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.213e-05
  |ΔD|_∞ = 2.669e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.1e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.708e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.426e-06
  |ΔD|_∞ = 5.339e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3468e-01 J
  → Fracture energy : 1.3585e+00 J
  → Total energy    : 1.4931e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0022.vtu


## Step 24/141: t = 5.91e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 23 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.15e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.360e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.875e-04
  |ΔD|_∞ = 7.158e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.15e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.175e-04
  |ΔD|_∞ = 1.432e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.15e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.350e-05
  |ΔD|_∞ = 2.863e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.15e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.700e-06
  |ΔD|_∞ = 5.726e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4711e-01 J
  → Fracture energy : 1.3585e+00 J
  → Total energy    : 1.5057e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0023.vtu


## Step 25/141: t = 6.17e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 24 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  ||Δu||/||u|| = 4.180e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.232e-04
  |ΔD|_∞ = 7.678e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.392e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.246e-04
  |ΔD|_∞ = 1.536e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.493e-05
  |ΔD|_∞ = 3.071e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.112e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.986e-06
  |ΔD|_∞ = 6.142e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6008e-01 J
  → Fracture energy : 1.3586e+00 J
  → Total energy    : 1.5187e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0024.vtu


## Step 26/141: t = 6.43e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 25 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.25e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.014e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.609e-04
  |ΔD|_∞ = 8.239e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.25e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.289e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.322e-04
  |ΔD|_∞ = 1.648e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.25e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.644e-05
  |ΔD|_∞ = 3.296e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.25e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.289e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.287e-06
  |ΔD|_∞ = 6.592e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7358e-01 J
  → Fracture energy : 1.3588e+00 J
  → Total energy    : 1.5323e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0025.vtu


## Step 27/141: t = 6.69e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 26 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.861e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.007e-04
  |ΔD|_∞ = 8.848e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.133e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.401e-04
  |ΔD|_∞ = 1.770e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.133e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.803e-05
  |ΔD|_∞ = 3.539e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.605e-06
  |ΔD|_∞ = 7.078e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8761e-01 J
  → Fracture energy : 1.3589e+00 J
  → Total energy    : 1.5465e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0026.vtu


## Step 28/141: t = 6.94e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 27 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.35e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.720e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.427e-04
  |ΔD|_∞ = 9.512e-04

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.35e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.931e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.485e-04
  |ΔD|_∞ = 1.902e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.35e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.931e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.971e-05
  |ΔD|_∞ = 3.805e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.35e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.942e-06
  |ΔD|_∞ = 7.609e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0216e-01 J
  → Fracture energy : 1.3590e+00 J
  → Total energy    : 1.5612e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0027.vtu


## Step 29/141: t = 7.20e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 28 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  ||Δu||/||u|| = 3.590e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.876e-04
  |ΔD|_∞ = 1.024e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.286e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.575e-04
  |ΔD|_∞ = 2.048e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.754e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.150e-05
  |ΔD|_∞ = 4.096e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.236e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.300e-06
  |ΔD|_∞ = 8.192e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1723e-01 J
  → Fracture energy : 1.3592e+00 J
  → Total energy    : 1.5764e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0028.vtu


## Step 30/141: t = 7.46e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 29 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.45e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.468e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.354e-04
  |ΔD|_∞ = 1.105e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.45e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.671e-04
  |ΔD|_∞ = 2.209e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.45e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.342e-05
  |ΔD|_∞ = 4.418e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.45e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.684e-06
  |ΔD|_∞ = 8.837e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3282e-01 J
  → Fracture energy : 1.3594e+00 J
  → Total energy    : 1.5922e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0029.vtu


## Step 31/141: t = 7.71e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 30 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.5e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.355e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.871e-04
  |ΔD|_∞ = 1.194e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.5e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.566e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.774e-04
  |ΔD|_∞ = 2.388e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.5e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.232e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.548e-05
  |ΔD|_∞ = 4.776e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.5e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.242e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.097e-06
  |ΔD|_∞ = 9.553e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4891e-01 J
  → Fracture energy : 1.3596e+00 J
  → Total energy    : 1.6085e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0030.vtu


## Step 32/141: t = 7.97e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 31 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.55e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.250e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.430e-04
  |ΔD|_∞ = 1.295e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.55e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.771e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.886e-04
  |ΔD|_∞ = 2.589e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.55e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.772e-05
  |ΔD|_∞ = 5.178e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.55e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.544e-06
  |ΔD|_∞ = 1.036e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6550e-01 J
  → Fracture energy : 1.3599e+00 J
  → Total energy    : 1.6254e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0031.vtu


## Step 33/141: t = 8.23e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 32 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  ||Δu||/||u|| = 3.152e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.004e-03
  |ΔD|_∞ = 1.408e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.608e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.008e-04
  |ΔD|_∞ = 2.816e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.739e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.016e-05
  |ΔD|_∞ = 5.632e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.901e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.032e-06
  |ΔD|_∞ = 1.126e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8258e-01 J
  → Fracture energy : 1.3602e+00 J
  → Total energy    : 1.6428e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0032.vtu


## Step 34/141: t = 8.49e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 33 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.65e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.061e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.071e-03
  |ΔD|_∞ = 1.538e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.65e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.238e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.143e-04
  |ΔD|_∞ = 3.076e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.65e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.238e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.285e-05
  |ΔD|_∞ = 6.152e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.65e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.570e-06
  |ΔD|_∞ = 1.230e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.0014e-01 J
  → Fracture energy : 1.3605e+00 J
  → Total energy    : 1.6607e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0033.vtu


## Step 35/141: t = 8.74e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 34 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.7e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.975e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.146e-03
  |ΔD|_∞ = 1.688e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.7e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.187e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.291e-04
  |ΔD|_∞ = 3.375e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.7e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.583e-05
  |ΔD|_∞ = 6.750e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.7e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.166e-06
  |ΔD|_∞ = 1.350e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1817e-01 J
  → Fracture energy : 1.3610e+00 J
  → Total energy    : 1.6791e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0034.vtu


## Step 36/141: t = 9.00e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 35 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.75e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.896e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.229e-03
  |ΔD|_∞ = 1.863e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.75e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.459e-04
  |ΔD|_∞ = 3.726e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.75e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.918e-05
  |ΔD|_∞ = 7.452e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.75e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.836e-06
  |ΔD|_∞ = 1.490e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3666e-01 J
  → Fracture energy : 1.3614e+00 J
  → Total energy    : 1.6981e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0035.vtu


## Step 37/141: t = 9.26e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 36 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  ||Δu||/||u|| = 2.822e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.325e-03
  |ΔD|_∞ = 2.070e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.278e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.650e-04
  |ΔD|_∞ = 4.140e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.278e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.299e-05
  |ΔD|_∞ = 8.281e-05

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
  ||ΔD||/||D|| = 1.060e-05
  |ΔD|_∞ = 1.656e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.120e-06
  |ΔD|_∞ = 3.312e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.5558e-01 J
  → Fracture energy : 1.3620e+00 J
  → Total energy    : 1.7176e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0036.vtu


## Step 38/141: t = 9.51e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 37 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.85e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.754e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.434e-03
  |ΔD|_∞ = 2.317e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.85e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.311e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.867e-04
  |ΔD|_∞ = 4.635e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.85e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.734e-05
  |ΔD|_∞ = 9.270e-05

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.85e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.147e-05
  |ΔD|_∞ = 1.854e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.85e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.294e-06
  |ΔD|_∞ = 3.708e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.7491e-01 J
  → Fracture energy : 1.3626e+00 J
  → Total energy    : 1.7375e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0037.vtu


## Step 39/141: t = 9.77e+02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 38 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.9e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.691e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.562e-03
  |ΔD|_∞ = 2.619e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.9e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.125e-04
  |ΔD|_∞ = 5.238e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.9e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.249e-05
  |ΔD|_∞ = 1.048e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.9e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.250e-05
  |ΔD|_∞ = 2.095e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.9e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.500e-06
  |ΔD|_∞ = 4.191e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.9462e-01 J
  → Fracture energy : 1.3634e+00 J
  → Total energy    : 1.7580e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0038.vtu


## Step 40/141: t = 1.00e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 39 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.95e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.635e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.717e-03
  |ΔD|_∞ = 2.995e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.95e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.516e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.434e-04
  |ΔD|_∞ = 5.989e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.95e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.868e-05
  |ΔD|_∞ = 1.198e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.95e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.047e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.374e-05
  |ΔD|_∞ = 2.396e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 1.95e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.047e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.747e-06
  |ΔD|_∞ = 4.791e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.1466e-01 J
  → Fracture energy : 1.3643e+00 J
  → Total energy    : 1.7790e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0039.vtu


## Step 41/141: t = 1.03e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 40 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  ||Δu||/||u|| = 2.586e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.908e-03
  |ΔD|_∞ = 3.467e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.813e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.815e-04
  |ΔD|_∞ = 6.934e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.973e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.630e-05
  |ΔD|_∞ = 1.387e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.972e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.526e-05
  |ΔD|_∞ = 2.773e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.052e-06
  |ΔD|_∞ = 5.547e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.3498e-01 J
  → Fracture energy : 1.3655e+00 J
  → Total energy    : 1.8004e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0040.vtu


## Step 42/141: t = 1.05e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 41 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.05e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.547e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.149e-03
  |ΔD|_∞ = 4.073e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.05e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.474e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.299e-04
  |ΔD|_∞ = 8.146e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.05e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.598e-05
  |ΔD|_∞ = 1.629e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.05e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.720e-05
  |ΔD|_∞ = 3.259e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.05e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.476e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.439e-06
  |ΔD|_∞ = 6.517e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.5550e-01 J
  → Fracture energy : 1.3668e+00 J
  → Total energy    : 1.8223e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0041.vtu


## Step 43/141: t = 1.08e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 42 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.1e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.520e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.468e-03
  |ΔD|_∞ = 4.864e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.1e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.100e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.937e-04
  |ΔD|_∞ = 9.728e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.1e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.248e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.873e-05
  |ΔD|_∞ = 1.946e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.1e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.325e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.975e-05
  |ΔD|_∞ = 3.891e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.1e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.645e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.949e-06
  |ΔD|_∞ = 7.783e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.7606e-01 J
  → Fracture energy : 1.3686e+00 J
  → Total energy    : 1.8447e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0042.vtu


## Step 44/141: t = 1.11e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 43 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.15e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.515e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.907e-03
  |ΔD|_∞ = 5.907e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.15e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.279e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.815e-04
  |ΔD|_∞ = 1.181e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.15e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.311e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.163e-04
  |ΔD|_∞ = 2.363e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.15e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.901e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.326e-05
  |ΔD|_∞ = 4.726e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.15e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.652e-06
  |ΔD|_∞ = 9.452e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.9646e-01 J
  → Fracture energy : 1.3709e+00 J
  → Total energy    : 1.8673e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0043.vtu


## Step 45/141: t = 1.13e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 44 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  ||Δu||/||u|| = 2.547e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.543e-03
  |ΔD|_∞ = 7.277e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.522e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.086e-04
  |ΔD|_∞ = 1.455e-03

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
  ||ΔD||/||D|| = 1.417e-04
  |ΔD|_∞ = 2.911e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.834e-05
  |ΔD|_∞ = 5.822e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.668e-06
  |ΔD|_∞ = 1.164e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.1628e-01 J
  → Fracture energy : 1.3740e+00 J
  → Total energy    : 1.8902e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0044.vtu


## Step 46/141: t = 1.16e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 45 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.25e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.658e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.507e-03
  |ΔD|_∞ = 9.657e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.25e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.013e-04
  |ΔD|_∞ = 1.931e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.25e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.803e-04
  |ΔD|_∞ = 3.863e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.25e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.553e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.605e-05
  |ΔD|_∞ = 7.726e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.25e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.553e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.211e-06
  |ΔD|_∞ = 1.545e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.3474e-01 J
  → Fracture energy : 1.3784e+00 J
  → Total energy    : 1.9131e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0045.vtu


## Step 47/141: t = 1.18e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 46 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.954e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.007e-03
  |ΔD|_∞ = 1.375e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.192e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.201e-03
  |ΔD|_∞ = 2.750e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.306e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.403e-04
  |ΔD|_∞ = 5.500e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.041e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.806e-05
  |ΔD|_∞ = 1.100e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.612e-06
  |ΔD|_∞ = 2.200e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.5019e-01 J
  → Fracture energy : 1.3853e+00 J
  → Total energy    : 1.9354e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0046.vtu


## Step 48/141: t = 1.21e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 47 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.35e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.729e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.247e-03
  |ΔD|_∞ = 2.015e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.35e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.292e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.649e-03
  |ΔD|_∞ = 4.029e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.35e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.116e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.299e-04
  |ΔD|_∞ = 8.058e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.35e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.598e-05
  |ΔD|_∞ = 1.612e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.35e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.152e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.320e-05
  |ΔD|_∞ = 3.223e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.35e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.153e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.639e-06
  |ΔD|_∞ = 6.447e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.5893e-01 J
  → Fracture energy : 1.3965e+00 J
  → Total energy    : 1.9554e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0047.vtu


## Step 49/141: t = 1.23e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 48 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  ||Δu||/||u|| = 5.745e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.098e-02
  |ΔD|_∞ = 2.592e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.215e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.196e-03
  |ΔD|_∞ = 5.185e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.392e-04
  |ΔD|_∞ = 1.037e-03

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
  ||ΔD||/||D|| = 8.784e-05
  |ΔD|_∞ = 2.074e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.168e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.757e-05
  |ΔD|_∞ = 4.148e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.838e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.514e-06
  |ΔD|_∞ = 8.296e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.5344e-01 J
  → Fracture energy : 1.4148e+00 J
  → Total energy    : 1.9682e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0048.vtu


## Step 50/141: t = 1.26e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 49 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.45e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.032e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.241e-02
  |ΔD|_∞ = 3.359e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.45e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.354e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.483e-03
  |ΔD|_∞ = 6.718e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.45e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.965e-04
  |ΔD|_∞ = 1.344e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.45e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.088e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.930e-05
  |ΔD|_∞ = 2.687e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.45e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.088e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.986e-05
  |ΔD|_∞ = 5.375e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.45e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.354e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.972e-06
  |ΔD|_∞ = 1.075e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.2455e-01 J
  → Fracture energy : 1.4396e+00 J
  → Total energy    : 1.9641e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0049.vtu


## Step 51/141: t = 1.29e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 50 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.5e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.611e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.712e-03
  |ΔD|_∞ = 3.764e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.5e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.942e-03
  |ΔD|_∞ = 7.529e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.5e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.884e-04
  |ΔD|_∞ = 1.506e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.5e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.107e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.769e-05
  |ΔD|_∞ = 3.012e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.5e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.107e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.554e-05
  |ΔD|_∞ = 6.023e-05

Convergence check


#### Iteration 6/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.5e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.107e-06
  |ΔD|_∞ = 1.205e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 6 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.8478e-01 J
  → Fracture energy : 1.4609e+00 J
  → Total energy    : 1.9457e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0050.vtu


## Step 52/141: t = 1.31e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 51 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.55e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.226e-01

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.464e-03
  |ΔD|_∞ = 2.279e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.55e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.093e-03
  |ΔD|_∞ = 4.558e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.55e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.811e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.186e-04
  |ΔD|_∞ = 9.116e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.55e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.811e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.371e-05
  |ΔD|_∞ = 1.823e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.55e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.086e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.742e-06
  |ΔD|_∞ = 3.646e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.7547e-01 J
  → Fracture energy : 1.4714e+00 J
  → Total energy    : 1.9468e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0051.vtu


## Step 53/141: t = 1.34e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 52 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  ||Δu||/||u|| = 4.449e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.328e-03
  |ΔD|_∞ = 1.955e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.761e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.655e-04
  |ΔD|_∞ = 3.910e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.761e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.331e-04
  |ΔD|_∞ = 7.820e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.662e-05
  |ΔD|_∞ = 1.564e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.732e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.324e-06
  |ΔD|_∞ = 3.128e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.8610e-01 J
  → Fracture energy : 1.4759e+00 J
  → Total energy    : 1.9620e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0052.vtu


## Step 54/141: t = 1.36e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 53 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.65e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.565e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.244e-03
  |ΔD|_∞ = 1.657e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.65e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.091e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.487e-04
  |ΔD|_∞ = 3.313e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.65e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.996e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.974e-05
  |ΔD|_∞ = 6.626e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.65e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.996e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.795e-05
  |ΔD|_∞ = 1.325e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.65e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.964e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.590e-06
  |ΔD|_∞ = 2.651e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.0150e-01 J
  → Fracture energy : 1.4784e+00 J
  → Total energy    : 1.9799e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0053.vtu


## Step 55/141: t = 1.39e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 54 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.7e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.120e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.605e-03
  |ΔD|_∞ = 1.187e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.7e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.210e-04
  |ΔD|_∞ = 2.373e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.7e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.421e-05
  |ΔD|_∞ = 4.747e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.7e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.853e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.284e-05
  |ΔD|_∞ = 9.494e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.7e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.035e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.568e-06
  |ΔD|_∞ = 1.899e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.1863e-01 J
  → Fracture energy : 1.4799e+00 J
  → Total energy    : 1.9986e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0054.vtu


## Step 56/141: t = 1.41e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 55 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.75e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.957e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.193e-03
  |ΔD|_∞ = 1.117e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.75e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.908e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.385e-04
  |ΔD|_∞ = 2.234e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.75e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.419e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.771e-05
  |ΔD|_∞ = 4.468e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.75e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.406e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.542e-06
  |ΔD|_∞ = 8.936e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.3673e-01 J
  → Fracture energy : 1.4811e+00 J
  → Total energy    : 2.0178e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0055.vtu


## Step 57/141: t = 1.44e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 56 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  ||Δu||/||u|| = 1.871e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.382e-04
  |ΔD|_∞ = 9.246e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.876e-04
  |ΔD|_∞ = 1.849e-03

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
  ||ΔD||/||D|| = 3.753e-05
  |ΔD|_∞ = 3.699e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.474e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.506e-06
  |ΔD|_∞ = 7.397e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.5546e-01 J
  → Fracture energy : 1.4820e+00 J
  → Total energy    : 2.0374e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0056.vtu


## Step 58/141: t = 1.47e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 57 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.85e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.814e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.917e-04
  |ΔD|_∞ = 7.111e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.85e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.637e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.583e-04
  |ΔD|_∞ = 1.422e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.85e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.167e-05
  |ΔD|_∞ = 2.844e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.85e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.978e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.334e-06
  |ΔD|_∞ = 5.689e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.7467e-01 J
  → Fracture energy : 1.4827e+00 J
  → Total energy    : 2.0574e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0057.vtu


## Step 59/141: t = 1.49e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 58 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.9e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.771e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.127e-04
  |ΔD|_∞ = 5.467e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.9e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.425e-04
  |ΔD|_∞ = 1.093e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.9e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.332e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.851e-05
  |ΔD|_∞ = 2.187e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.9e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.332e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.702e-06
  |ΔD|_∞ = 4.374e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.9430e-01 J
  → Fracture energy : 1.4834e+00 J
  → Total energy    : 2.0777e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0058.vtu


## Step 60/141: t = 1.52e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 59 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.95e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.735e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.722e-04
  |ΔD|_∞ = 5.116e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.95e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.988e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.344e-04
  |ΔD|_∞ = 1.023e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.95e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.088e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.689e-05
  |ΔD|_∞ = 2.047e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 2.95e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.351e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.378e-06
  |ΔD|_∞ = 4.093e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.1431e-01 J
  → Fracture energy : 1.4841e+00 J
  → Total energy    : 2.0984e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0059.vtu


## Step 61/141: t = 1.54e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 60 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  ||Δu||/||u|| = 1.702e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.518e-04
  |ΔD|_∞ = 4.750e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.382e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.304e-04
  |ΔD|_∞ = 9.501e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.607e-05
  |ΔD|_∞ = 1.900e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.214e-06
  |ΔD|_∞ = 3.800e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.3467e-01 J
  → Fracture energy : 1.4847e+00 J
  → Total energy    : 2.1193e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0060.vtu


## Step 62/141: t = 1.57e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 61 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.05e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.671e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.416e-04
  |ΔD|_∞ = 4.339e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.05e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.645e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.283e-04
  |ΔD|_∞ = 8.678e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.05e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.645e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.566e-05
  |ΔD|_∞ = 1.736e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.05e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.133e-06
  |ΔD|_∞ = 3.471e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.5539e-01 J
  → Fracture energy : 1.4853e+00 J
  → Total energy    : 2.1407e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0061.vtu


## Step 63/141: t = 1.59e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 62 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.1e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.642e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.393e-04
  |ΔD|_∞ = 4.103e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.1e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.289e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.279e-04
  |ΔD|_∞ = 8.206e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.1e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.258e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.557e-05
  |ΔD|_∞ = 1.641e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.1e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.243e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.115e-06
  |ΔD|_∞ = 3.282e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.7645e-01 J
  → Fracture energy : 1.4859e+00 J
  → Total energy    : 2.1623e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0062.vtu


## Step 64/141: t = 1.62e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 63 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.15e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.614e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.393e-04
  |ΔD|_∞ = 4.197e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.15e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.431e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.279e-04
  |ΔD|_∞ = 8.394e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.15e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.557e-05
  |ΔD|_∞ = 1.679e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.15e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.114e-06
  |ΔD|_∞ = 3.358e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.9785e-01 J
  → Fracture energy : 1.4864e+00 J
  → Total energy    : 2.1843e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0063.vtu


## Step 65/141: t = 1.65e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 64 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  ||Δu||/||u|| = 1.587e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.411e-04
  |ΔD|_∞ = 4.209e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.431e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.282e-04
  |ΔD|_∞ = 8.417e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.481e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.564e-05
  |ΔD|_∞ = 1.683e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.129e-06
  |ΔD|_∞ = 3.367e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.1959e-01 J
  → Fracture energy : 1.4870e+00 J
  → Total energy    : 2.2066e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0064.vtu


## Step 66/141: t = 1.67e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 65 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.25e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.562e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.432e-04
  |ΔD|_∞ = 4.116e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.25e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.936e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.286e-04
  |ΔD|_∞ = 8.233e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.25e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.116e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.573e-05
  |ΔD|_∞ = 1.647e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.25e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.146e-06
  |ΔD|_∞ = 3.293e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.4166e-01 J
  → Fracture energy : 1.4876e+00 J
  → Total energy    : 2.2292e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0065.vtu


## Step 67/141: t = 1.70e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 66 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.537e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.446e-04
  |ΔD|_∞ = 3.964e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.158e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.289e-04
  |ΔD|_∞ = 7.927e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.578e-05
  |ΔD|_∞ = 1.585e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.157e-06
  |ΔD|_∞ = 3.171e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.6408e-01 J
  → Fracture energy : 1.4881e+00 J
  → Total energy    : 2.2522e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0066.vtu


## Step 68/141: t = 1.72e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 67 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.35e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.513e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.450e-04
  |ΔD|_∞ = 3.903e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.35e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.290e-04
  |ΔD|_∞ = 7.806e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.35e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.431e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.580e-05
  |ΔD|_∞ = 1.561e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.35e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.431e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.160e-06
  |ΔD|_∞ = 3.122e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.8683e-01 J
  → Fracture energy : 1.4887e+00 J
  → Total energy    : 2.2755e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0067.vtu


## Step 69/141: t = 1.75e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 68 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  ||Δu||/||u|| = 1.490e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.449e-04
  |ΔD|_∞ = 4.005e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.290e-04
  |ΔD|_∞ = 8.010e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.580e-05
  |ΔD|_∞ = 1.602e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.598e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.159e-06
  |ΔD|_∞ = 3.204e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.0992e-01 J
  → Fracture energy : 1.4893e+00 J
  → Total energy    : 2.2992e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0068.vtu


## Step 70/141: t = 1.77e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 69 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.45e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.468e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.449e-04
  |ΔD|_∞ = 4.104e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.45e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.167e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.290e-04
  |ΔD|_∞ = 8.208e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.45e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.167e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.579e-05
  |ΔD|_∞ = 1.642e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.45e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.682e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.159e-06
  |ΔD|_∞ = 3.283e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.3334e-01 J
  → Fracture energy : 1.4898e+00 J
  → Total energy    : 2.3231e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0069.vtu


## Step 71/141: t = 1.80e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 70 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.5e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.446e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.459e-04
  |ΔD|_∞ = 4.202e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.5e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.700e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.292e-04
  |ΔD|_∞ = 8.404e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.5e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.584e-05
  |ΔD|_∞ = 1.681e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.5e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.167e-06
  |ΔD|_∞ = 3.362e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.5710e-01 J
  → Fracture energy : 1.4904e+00 J
  → Total energy    : 2.3475e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0070.vtu


## Step 72/141: t = 1.83e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 71 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.55e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.425e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.480e-04
  |ΔD|_∞ = 4.302e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.55e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.362e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.296e-04
  |ΔD|_∞ = 8.603e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.55e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.592e-05
  |ΔD|_∞ = 1.721e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.55e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.184e-06
  |ΔD|_∞ = 3.441e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.8119e-01 J
  → Fracture energy : 1.4909e+00 J
  → Total energy    : 2.3721e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0071.vtu


## Step 73/141: t = 1.85e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 72 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  ||Δu||/||u|| = 1.404e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.516e-04
  |ΔD|_∞ = 4.405e-03

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
  ||ΔD||/||D|| = 1.303e-04
  |ΔD|_∞ = 8.810e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.606e-05
  |ΔD|_∞ = 1.762e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.213e-06
  |ΔD|_∞ = 3.524e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.0562e-01 J
  → Fracture energy : 1.4915e+00 J
  → Total energy    : 2.3971e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0072.vtu


## Step 74/141: t = 1.88e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 73 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.65e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.384e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.575e-04
  |ΔD|_∞ = 4.520e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.65e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.590e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.315e-04
  |ΔD|_∞ = 9.039e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.65e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.630e-05
  |ΔD|_∞ = 1.808e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.65e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.750e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.260e-06
  |ΔD|_∞ = 3.616e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.3038e-01 J
  → Fracture energy : 1.4920e+00 J
  → Total energy    : 2.4224e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0073.vtu


## Step 75/141: t = 1.90e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 74 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.7e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.365e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.655e-04
  |ΔD|_∞ = 4.643e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.7e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.484e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.331e-04
  |ΔD|_∞ = 9.286e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.7e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.267e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.662e-05
  |ΔD|_∞ = 1.857e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.7e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.581e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.324e-06
  |ΔD|_∞ = 3.714e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.5546e-01 J
  → Fracture energy : 1.4926e+00 J
  → Total energy    : 2.4480e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0074.vtu


## Step 76/141: t = 1.93e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 75 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.75e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.347e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.747e-04
  |ΔD|_∞ = 4.779e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.75e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.600e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.349e-04
  |ΔD|_∞ = 9.558e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.75e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.600e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.699e-05
  |ΔD|_∞ = 1.912e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.75e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.397e-06
  |ΔD|_∞ = 3.823e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.8087e-01 J
  → Fracture energy : 1.4931e+00 J
  → Total energy    : 2.4740e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0075.vtu


## Step 77/141: t = 1.95e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 76 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  ||Δu||/||u|| = 1.329e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.856e-04
  |ΔD|_∞ = 4.928e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.371e-04
  |ΔD|_∞ = 9.855e-04

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.742e-05
  |ΔD|_∞ = 1.971e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.713e-19

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.485e-06
  |ΔD|_∞ = 3.942e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0066e+00 J
  → Fracture energy : 1.4937e+00 J
  → Total energy    : 2.5003e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0076.vtu


## Step 78/141: t = 1.98e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 77 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.85e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.311e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.975e-04
  |ΔD|_∞ = 5.086e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.85e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.326e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.395e-04
  |ΔD|_∞ = 1.017e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.85e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.326e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.790e-05
  |ΔD|_∞ = 2.034e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.85e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.580e-06
  |ΔD|_∞ = 4.069e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0327e+00 J
  → Fracture energy : 1.4943e+00 J
  → Total energy    : 2.5269e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0077.vtu


## Step 79/141: t = 2.01e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 78 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.9e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.294e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.110e-04
  |ΔD|_∞ = 5.256e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.9e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.422e-04
  |ΔD|_∞ = 1.051e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.9e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.844e-05
  |ΔD|_∞ = 2.102e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.9e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.688e-06
  |ΔD|_∞ = 4.205e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0590e+00 J
  → Fracture energy : 1.4949e+00 J
  → Total energy    : 2.5539e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0078.vtu


## Step 80/141: t = 2.03e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 79 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.95e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.277e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.257e-04
  |ΔD|_∞ = 5.441e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.95e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.451e-04
  |ΔD|_∞ = 1.088e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.95e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.903e-05
  |ΔD|_∞ = 2.177e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 3.95e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.360e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.805e-06
  |ΔD|_∞ = 4.353e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0857e+00 J
  → Fracture energy : 1.4954e+00 J
  → Total energy    : 2.5812e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0079.vtu


## Step 81/141: t = 2.06e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 80 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  ||Δu||/||u|| = 1.261e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.414e-04
  |ΔD|_∞ = 5.642e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.483e-04
  |ΔD|_∞ = 1.128e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.966e-05
  |ΔD|_∞ = 2.257e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.037e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.931e-06
  |ΔD|_∞ = 4.513e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1128e+00 J
  → Fracture energy : 1.4961e+00 J
  → Total energy    : 2.6088e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0080.vtu


## Step 82/141: t = 2.08e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 81 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.05e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.246e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.581e-04
  |ΔD|_∞ = 5.854e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.05e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.950e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.516e-04
  |ΔD|_∞ = 1.171e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.05e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.032e-05
  |ΔD|_∞ = 2.342e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.05e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.065e-06
  |ΔD|_∞ = 4.684e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1401e+00 J
  → Fracture energy : 1.4967e+00 J
  → Total energy    : 2.6368e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0081.vtu


## Step 83/141: t = 2.11e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 82 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.1e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.230e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.752e-04
  |ΔD|_∞ = 6.080e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.1e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.251e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.550e-04
  |ΔD|_∞ = 1.216e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.1e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.101e-05
  |ΔD|_∞ = 2.432e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.1e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.201e-06
  |ΔD|_∞ = 4.864e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1677e+00 J
  → Fracture energy : 1.4973e+00 J
  → Total energy    : 2.6651e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0082.vtu


## Step 84/141: t = 2.13e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 83 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.15e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.216e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.935e-04
  |ΔD|_∞ = 6.319e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.15e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.134e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.587e-04
  |ΔD|_∞ = 1.264e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.15e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.174e-05
  |ΔD|_∞ = 2.528e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.15e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.133e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.348e-06
  |ΔD|_∞ = 5.055e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1957e+00 J
  → Fracture energy : 1.4980e+00 J
  → Total energy    : 2.6937e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0083.vtu


## Step 85/141: t = 2.16e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 84 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  ||Δu||/||u|| = 1.201e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.121e-04
  |ΔD|_∞ = 6.572e-03

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
  ||ΔD||/||D|| = 1.624e-04
  |ΔD|_∞ = 1.314e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.180e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.248e-05
  |ΔD|_∞ = 2.629e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.180e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.497e-06
  |ΔD|_∞ = 5.258e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2240e+00 J
  → Fracture energy : 1.4987e+00 J
  → Total energy    : 2.7226e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0084.vtu


## Step 86/141: t = 2.19e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 85 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.25e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.187e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.321e-04
  |ΔD|_∞ = 6.840e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.25e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.757e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.664e-04
  |ΔD|_∞ = 1.368e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.25e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.499e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.328e-05
  |ΔD|_∞ = 2.736e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.25e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.702e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.657e-06
  |ΔD|_∞ = 5.472e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2525e+00 J
  → Fracture energy : 1.4994e+00 J
  → Total energy    : 2.7519e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0085.vtu


## Step 87/141: t = 2.21e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 86 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.173e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.528e-04
  |ΔD|_∞ = 7.125e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.079e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.706e-04
  |ΔD|_∞ = 1.425e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.673e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.411e-05
  |ΔD|_∞ = 2.850e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.081e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.823e-06
  |ΔD|_∞ = 5.700e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2814e+00 J
  → Fracture energy : 1.5001e+00 J
  → Total energy    : 2.7815e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0086.vtu


## Step 88/141: t = 2.24e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 87 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.35e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.160e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.748e-04
  |ΔD|_∞ = 7.428e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.35e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.042e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.750e-04
  |ΔD|_∞ = 1.486e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.35e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.031e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.499e-05
  |ΔD|_∞ = 2.971e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.35e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.611e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.998e-06
  |ΔD|_∞ = 5.943e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3106e+00 J
  → Fracture energy : 1.5008e+00 J
  → Total energy    : 2.8114e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0087.vtu


## Step 89/141: t = 2.26e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 88 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  ||Δu||/||u|| = 1.147e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.981e-04
  |ΔD|_∞ = 7.754e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.796e-04
  |ΔD|_∞ = 1.551e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.863e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.592e-05
  |ΔD|_∞ = 3.102e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.393e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.185e-06
  |ΔD|_∞ = 6.203e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3401e+00 J
  → Fracture energy : 1.5016e+00 J
  → Total energy    : 2.8417e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0088.vtu


## Step 90/141: t = 2.29e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 89 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.45e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.134e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.233e-04
  |ΔD|_∞ = 8.107e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.45e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.794e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.847e-04
  |ΔD|_∞ = 1.621e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.45e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.794e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.693e-05
  |ΔD|_∞ = 3.243e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.45e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.387e-06
  |ΔD|_∞ = 6.486e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3699e+00 J
  → Fracture energy : 1.5023e+00 J
  → Total energy    : 2.8722e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0089.vtu


## Step 91/141: t = 2.31e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 90 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.5e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.121e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.507e-04
  |ΔD|_∞ = 8.491e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.5e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.305e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.901e-04
  |ΔD|_∞ = 1.698e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.5e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.305e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.803e-05
  |ΔD|_∞ = 3.396e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.5e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.605e-06
  |ΔD|_∞ = 6.793e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4000e+00 J
  → Fracture energy : 1.5031e+00 J
  → Total energy    : 2.9031e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0090.vtu


## Step 92/141: t = 2.34e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 91 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.55e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.109e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.806e-04
  |ΔD|_∞ = 8.911e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.55e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.961e-04
  |ΔD|_∞ = 1.782e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.55e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.156e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.922e-05
  |ΔD|_∞ = 3.564e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.55e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.091e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.845e-06
  |ΔD|_∞ = 7.129e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4304e+00 J
  → Fracture energy : 1.5040e+00 J
  → Total energy    : 2.9344e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0091.vtu


## Step 93/141: t = 2.37e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 92 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  ||Δu||/||u|| = 1.097e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.014e-03
  |ΔD|_∞ = 9.375e-03

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
  ||ΔD||/||D|| = 2.027e-04
  |ΔD|_∞ = 1.875e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.728e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.055e-05
  |ΔD|_∞ = 3.750e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.024e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.109e-06
  |ΔD|_∞ = 7.500e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4611e+00 J
  → Fracture energy : 1.5048e+00 J
  → Total energy    : 2.9659e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0092.vtu


## Step 94/141: t = 2.39e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 93 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.65e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.086e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.050e-03
  |ΔD|_∞ = 9.889e-03

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.65e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.101e-04
  |ΔD|_∞ = 1.978e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.65e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.850e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.202e-05
  |ΔD|_∞ = 3.956e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.65e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.850e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.403e-06
  |ΔD|_∞ = 7.911e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4921e+00 J
  → Fracture energy : 1.5057e+00 J
  → Total energy    : 2.9978e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0093.vtu


## Step 95/141: t = 2.42e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 94 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.7e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.074e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.092e-03
  |ΔD|_∞ = 1.046e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.7e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.959e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.183e-04
  |ΔD|_∞ = 2.093e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.7e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.959e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.366e-05
  |ΔD|_∞ = 4.186e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.7e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.733e-06
  |ΔD|_∞ = 8.371e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5233e+00 J
  → Fracture energy : 1.5067e+00 J
  → Total energy    : 3.0300e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0094.vtu


## Step 96/141: t = 2.44e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 95 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.75e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.063e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.138e-03
  |ΔD|_∞ = 1.111e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.75e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.660e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.276e-04
  |ΔD|_∞ = 2.222e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.75e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.552e-05
  |ΔD|_∞ = 4.444e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.75e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.035e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.105e-06
  |ΔD|_∞ = 8.888e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5549e+00 J
  → Fracture energy : 1.5076e+00 J
  → Total energy    : 3.0625e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0095.vtu


## Step 97/141: t = 2.47e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 96 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  ||Δu||/||u|| = 1.053e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.191e-03
  |ΔD|_∞ = 1.184e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.027e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.382e-04
  |ΔD|_∞ = 2.368e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.765e-05
  |ΔD|_∞ = 4.736e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.348e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.529e-06
  |ΔD|_∞ = 9.473e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5867e+00 J
  → Fracture energy : 1.5087e+00 J
  → Total energy    : 3.0953e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0096.vtu


## Step 98/141: t = 2.49e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 97 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.85e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.042e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.252e-03
  |ΔD|_∞ = 1.267e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.85e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.019e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.504e-04
  |ΔD|_∞ = 2.534e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.85e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.008e-05
  |ΔD|_∞ = 5.068e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.85e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.002e-05
  |ΔD|_∞ = 1.014e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.85e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.003e-06
  |ΔD|_∞ = 2.027e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6187e+00 J
  → Fracture energy : 1.5097e+00 J
  → Total energy    : 3.1285e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0097.vtu


## Step 99/141: t = 2.52e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 98 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.9e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.032e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.321e-03
  |ΔD|_∞ = 1.361e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.9e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.241e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.642e-04
  |ΔD|_∞ = 2.722e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.9e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.241e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.284e-05
  |ΔD|_∞ = 5.444e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.9e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.415e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.057e-05
  |ΔD|_∞ = 1.089e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.9e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.715e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.114e-06
  |ΔD|_∞ = 2.178e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6511e+00 J
  → Fracture energy : 1.5109e+00 J
  → Total energy    : 3.1619e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0098.vtu


## Step 100/141: t = 2.55e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 99 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.95e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.023e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.400e-03
  |ΔD|_∞ = 1.468e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.95e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.371e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.801e-04
  |ΔD|_∞ = 2.935e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.95e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 7.487e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.602e-05
  |ΔD|_∞ = 5.871e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.95e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.201e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.120e-05
  |ΔD|_∞ = 1.174e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 4.95e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.506e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.241e-06
  |ΔD|_∞ = 2.348e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6836e+00 J
  → Fracture energy : 1.5120e+00 J
  → Total energy    : 3.1957e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0099.vtu


## Step 101/141: t = 2.57e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 100 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  ||Δu||/||u|| = 1.013e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.493e-03
  |ΔD|_∞ = 1.590e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.800e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.986e-04
  |ΔD|_∞ = 3.179e-03

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
  ||ΔD||/||D|| = 5.973e-05
  |ΔD|_∞ = 6.359e-04

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
  ||ΔD||/||D|| = 1.195e-05
  |ΔD|_∞ = 1.272e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.389e-06
  |ΔD|_∞ = 2.544e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7165e+00 J
  → Fracture energy : 1.5133e+00 J
  → Total energy    : 3.2298e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0100.vtu


## Step 102/141: t = 2.60e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 101 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.05e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.004e-02

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.603e-03
  |ΔD|_∞ = 1.730e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.05e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.220e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.206e-04
  |ΔD|_∞ = 3.460e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.05e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.252e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.412e-05
  |ΔD|_∞ = 6.920e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.05e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.661e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.282e-05
  |ΔD|_∞ = 1.384e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.05e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.670e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.565e-06
  |ΔD|_∞ = 2.768e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7495e+00 J
  → Fracture energy : 1.5147e+00 J
  → Total energy    : 3.2642e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0101.vtu


## Step 103/141: t = 2.62e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 102 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.1e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.958e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.732e-03
  |ΔD|_∞ = 1.894e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.1e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.483e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.465e-04
  |ΔD|_∞ = 3.789e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.1e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.741e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.929e-05
  |ΔD|_∞ = 7.577e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.1e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.741e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.386e-05
  |ΔD|_∞ = 1.515e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.1e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.772e-06
  |ΔD|_∞ = 3.031e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7827e+00 J
  → Fracture energy : 1.5162e+00 J
  → Total energy    : 3.2988e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0102.vtu


## Step 104/141: t = 2.65e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 103 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.15e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.877e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.880e-03
  |ΔD|_∞ = 2.098e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.15e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.761e-04
  |ΔD|_∞ = 4.197e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.15e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.189e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 7.522e-05
  |ΔD|_∞ = 8.393e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.15e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.504e-05
  |ΔD|_∞ = 1.679e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.15e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.009e-06
  |ΔD|_∞ = 3.357e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8161e+00 J
  → Fracture energy : 1.5178e+00 J
  → Total energy    : 3.3338e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0103.vtu


## Step 105/141: t = 2.67e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 104 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  ||Δu||/||u|| = 9.802e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.041e-03
  |ΔD|_∞ = 2.323e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.082e-04
  |ΔD|_∞ = 4.646e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.163e-05
  |ΔD|_∞ = 9.292e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.842e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.633e-05
  |ΔD|_∞ = 1.858e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.2e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 6.842e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.265e-06
  |ΔD|_∞ = 3.717e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8496e+00 J
  → Fracture energy : 1.5195e+00 J
  → Total energy    : 3.3691e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0104.vtu


## Step 106/141: t = 2.70e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 105 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.25e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.731e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.210e-03
  |ΔD|_∞ = 2.557e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.25e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.421e-04
  |ΔD|_∞ = 5.115e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.25e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 8.841e-05
  |ΔD|_∞ = 1.023e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.25e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.768e-05
  |ΔD|_∞ = 2.046e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.25e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.537e-06
  |ΔD|_∞ = 4.092e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8833e+00 J
  → Fracture energy : 1.5214e+00 J
  → Total energy    : 3.4046e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0105.vtu


## Step 107/141: t = 2.73e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 106 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.664e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.382e-03
  |ΔD|_∞ = 2.782e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.365e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.763e-04
  |ΔD|_∞ = 5.565e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.283e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 9.526e-05
  |ΔD|_∞ = 1.113e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.086e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.905e-05
  |ΔD|_∞ = 2.226e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.3e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.282e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.811e-06
  |ΔD|_∞ = 4.452e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9170e+00 J
  → Fracture energy : 1.5234e+00 J
  → Total energy    : 3.4404e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0106.vtu


## Step 108/141: t = 2.75e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 107 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.35e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.600e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.546e-03
  |ΔD|_∞ = 3.026e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.35e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.134e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.092e-04
  |ΔD|_∞ = 6.052e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.35e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.134e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.018e-04
  |ΔD|_∞ = 1.210e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.35e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.226e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.037e-05
  |ΔD|_∞ = 2.421e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.35e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.226e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.074e-06
  |ΔD|_∞ = 4.841e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9509e+00 J
  → Fracture energy : 1.5256e+00 J
  → Total energy    : 3.4765e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0107.vtu


## Step 109/141: t = 2.78e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 108 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  ||Δu||/||u|| = 9.536e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.706e-03
  |ΔD|_∞ = 3.280e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.039e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.411e-04
  |ΔD|_∞ = 6.560e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.082e-04
  |ΔD|_∞ = 1.312e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.039e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.164e-05
  |ΔD|_∞ = 2.624e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.4e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.421e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.329e-06
  |ΔD|_∞ = 5.248e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9849e+00 J
  → Fracture energy : 1.5280e+00 J
  → Total energy    : 3.5128e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0108.vtu


## Step 110/141: t = 2.80e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 109 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.45e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.475e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.850e-03
  |ΔD|_∞ = 3.474e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.45e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.701e-04
  |ΔD|_∞ = 6.948e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.45e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.140e-04
  |ΔD|_∞ = 1.390e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.45e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.280e-05
  |ΔD|_∞ = 2.779e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.45e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.560e-06
  |ΔD|_∞ = 5.558e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0189e+00 J
  → Fracture energy : 1.5305e+00 J
  → Total energy    : 3.5494e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0109.vtu


## Step 111/141: t = 2.83e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 110 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.5e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.415e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.978e-03
  |ΔD|_∞ = 3.733e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.5e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.955e-04
  |ΔD|_∞ = 7.466e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.5e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.797e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.191e-04
  |ΔD|_∞ = 1.493e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.5e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.337e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.382e-05
  |ΔD|_∞ = 2.986e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.5e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 3.424e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.764e-06
  |ΔD|_∞ = 5.973e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0531e+00 J
  → Fracture energy : 1.5331e+00 J
  → Total energy    : 3.5862e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0110.vtu


## Step 112/141: t = 2.85e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 111 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.55e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.355e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.103e-03
  |ΔD|_∞ = 3.949e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.55e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.827e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.207e-04
  |ΔD|_∞ = 7.897e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.55e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.241e-04
  |ΔD|_∞ = 1.579e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.55e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.726e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.483e-05
  |ΔD|_∞ = 3.159e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.55e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.726e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 4.965e-06
  |ΔD|_∞ = 6.318e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0873e+00 J
  → Fracture energy : 1.5359e+00 J
  → Total energy    : 3.6232e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0111.vtu


## Step 113/141: t = 2.88e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 112 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  ||Δu||/||u|| = 9.299e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.207e-03
  |ΔD|_∞ = 4.143e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.414e-04
  |ΔD|_∞ = 8.286e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.283e-04
  |ΔD|_∞ = 1.657e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.566e-05
  |ΔD|_∞ = 3.314e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.6e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.131e-06
  |ΔD|_∞ = 6.629e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1216e+00 J
  → Fracture energy : 1.5389e+00 J
  → Total energy    : 3.6605e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0112.vtu


## Step 114/141: t = 2.91e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 113 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.65e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.241e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.310e-03
  |ΔD|_∞ = 4.383e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.65e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 4.691e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.620e-04
  |ΔD|_∞ = 8.767e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.65e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.544e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.324e-04
  |ΔD|_∞ = 1.753e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.65e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 2.517e-16

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.648e-05
  |ΔD|_∞ = 3.507e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.65e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.296e-06
  |ΔD|_∞ = 7.013e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1559e+00 J
  → Fracture energy : 1.5421e+00 J
  → Total energy    : 3.6980e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0113.vtu


## Step 115/141: t = 2.93e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 114 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.7e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.185e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.395e-03
  |ΔD|_∞ = 4.472e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.7e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.140e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.791e-04
  |ΔD|_∞ = 8.943e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.7e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.140e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.358e-04
  |ΔD|_∞ = 1.789e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.7e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.716e-05
  |ΔD|_∞ = 3.577e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.7e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 8.123e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.433e-06
  |ΔD|_∞ = 7.155e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1903e+00 J
  → Fracture energy : 1.5454e+00 J
  → Total energy    : 3.7357e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0114.vtu


## Step 116/141: t = 2.96e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 115 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-05
  → Staggering tolerance |ΔD|   : 1.0e-05
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-05
  → Relative tolerance dmg      : 1.0e-05


#### Iteration 1/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.75e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 9.128e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.480e-03
  |ΔD|_∞ = 4.690e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.75e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 5.544e-18

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 6.961e-04
  |ΔD|_∞ = 9.379e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.75e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 1.392e-04
  |ΔD|_∞ = 1.876e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.75e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 0.000e+00

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 2.784e-05
  |ΔD|_∞ = 3.752e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.75e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
  Linear solver
  ||Δu||/||u|| = 1.433e-17

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 5.569e-06
  |ΔD|_∞ = 7.504e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2248e+00 J
  → Fracture energy : 1.5488e+00 J
  → Total energy    : 3.7736e+00 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0115.vtu


## Step 117/141: t = 2.98e+03 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 116 | dt = 2.57e+01 s
Coupling = staggered
  → Max iterations              : 200
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
  ||Δu||/||u|| = 9.072e-03

**[INFO]** Assembling damage (AT2) problem...
Solving damage problem for 'steel' material
  ||ΔD||/||D|| = 3.559e-03
  |ΔD|_∞ = 4.758e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 1 → [0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 3 → [0.0, 5.8e-06]
  Building weak form, volume integrals (dx) for steel, tag = 6
