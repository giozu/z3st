Info    : Reading 'mesh.msh'...
Info    : 27 entities
Info    : 125 nodes
Info    : 160 elements
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
  → Time steps          : 41
  → Regime              : 3d
  → Models active       :
      thermal    → OFF
      mechanical → ON
      damage     → OFF
      cluster    → OFF
      plasticity → ON
  → Gap conductance     : None (value = 0.0)



### FiniteElementSetup initializer

Mechanical element order: 1
Thermal function space (V_t): FunctionSpace(<Mesh #0>, Basix element (P, hexahedron, 1, gll_warped, unset, False, float64, []))
Mechanical function space (V_m): FunctionSpace(<Mesh #0>, blocked element (Basix element (P, hexahedron, 1, gll_warped, unset, False, float64, []), (3,)))
Scalar function space (Q): FunctionSpace(<Mesh #0>, Basix element (P, hexahedron, 0, gll_warped, unset, True, float64, []))
Plasticity function spaces (V_pl_tensor, Q_pl) initializing with Quadrature degree 2
Plasticity function space (V_pl_tensor): FunctionSpace(<Mesh #0>, blocked element (QuadratureElement(hexahedron, array([[0.21132487, 0.21132487, 0.21132487],       [0.21132487, 0.21132487, 0.78867513],       [0.21132487, 0.78867513, 0.21132487],       [0.21132487, 0.78867513, 0.78867513],       [0.78867513, 0.21132487, 0.21132487],       [0.78867513, 0.21132487, 0.78867513],       [0.78867513, 0.78867513, 0.21132487],       [0.78867513, 0.78867513, 0.78867513]]), array([0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125]), IdentityPullback()), (3, 3)))
Plasticity function space (Q_pl): FunctionSpace(<Mesh #0>, QuadratureElement(hexahedron, array([[0.21132487, 0.21132487, 0.21132487],       [0.21132487, 0.21132487, 0.78867513],       [0.21132487, 0.78867513, 0.21132487],       [0.21132487, 0.78867513, 0.78867513],       [0.78867513, 0.21132487, 0.21132487],       [0.78867513, 0.21132487, 0.78867513],       [0.78867513, 0.78867513, 0.21132487],       [0.78867513, 0.78867513, 0.78867513]]), array([0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125]), IdentityPullback()))
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
  rtol                : 1e-07
  convergence         : rel_norm
[PlasticityModel] initializer
[spine.load_materials]
Material loaded: grain
  → k defined as constant: 50.0
  → Gc not defined for grain
  → constitutive model: custom
  E               → 200000000000.0 (float)
  G               → 76923076923.07692 (float)
  T_ref           → 300.0 (float)
  alpha           → 1.2e-05 (float)
  bulk_modulus    → 166666666666.66666 (float)
  constitutive    → custom (str)
  constitutive_mode → custom (str)
  cp              → 400.0 (float)
  g0              → 200000000.0 (float)
  gamma0          → 0.001 (float)
  k               → 50.0 (float)
  lmbda           → 115384615384.61539 (float)
  n_pow           → 5.0 (float)
  name            → single_crystal (str)
  nu              → 0.3 (float)
  rho             → 7800.0 (float)
  stress_function → single_crystal_law.single_crystal_stress (str)
[spine.initialize_fields]

Initializing the displacement field...
  Initial u: min=0.00e+00 m, max=0.00e+00 m, mean=0.00e+00 m



***


### spine - set_boundary_conditions


***



Loading boundary conditions from 'boundary_conditions.yaml'
  **[INFO]** Constant Dirichlet vector (3D) → [0.0, 0.0, 0.0]
  **[INFO]** Dirichlet mechanical BC on 'grain' → [0.0, 0.0, 0.0] at region 'zmin'
  **[INFO]** Dirichlet_z mechanical BC on 'grain' → 0.0 (first step) at region 'zmax'


## Step 01/41: t = 0.00e+00 s | LHR = 0.00e+00 W/m

  → dt=0: solving static step / initial condition
Computing symbolic result fields (strain, stress, ...)

======================================================================
CRYSTAL PLASTICITY - SCHMID FACTOR CALCULATION
======================================================================
Slip system: (111)[0-11]
  Plane normal n = [0.57735027 0.57735027 0.57735027]
  Slip direction m = [ 0.          0.70710678 -0.70710678]

Loading direction: e_z = [0, 0, 1]

Schmid factor calculations:
  Method 1 (direct):        μ = |m·e_z| × |n·e_z| = 0.408248
  Method 2 (tensor P_zz):   μ = |P_zz|           = 0.408248
  Method 3 (full tensor):   μ = |P:σ_zz|         = 0.408248

Schmid tensor P:
  [ 0.          0.20412415 -0.20412415]
  [0.20412415 0.40824829 0.        ]
  [-0.20412415  0.         -0.40824829]

For uniaxial stress σ_zz:
  Resolved shear stress: τ = 0.408248 × σ_zz
======================================================================




***


### spine - solve


***



Current step = 0 | dt = 0.00e+00 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.0
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 1 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0000.vtu


## Step 02/41: t = 5.00e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 1 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00025
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.000e+00

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00025
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.521e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0001.vtu


## Step 03/41: t = 1.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 2 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.0005
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.000e-01

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.0005
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.120e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0002.vtu


## Step 04/41: t = 1.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 3 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00075
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.333e-01

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00075
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.167e-09

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0003.vtu


## Step 05/41: t = 2.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 4 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.001
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.500e-01

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.001
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 7.212e-09

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0004.vtu


## Step 06/41: t = 2.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 5 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00125
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.000e-01

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00125
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.198e-08

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0005.vtu


## Step 07/41: t = 3.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 6 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.0015
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.667e-01

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.0015
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.134e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0006.vtu


## Step 08/41: t = 3.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 7 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00175
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.428e-01

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00175
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.265e-14

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0007.vtu


## Step 09/41: t = 4.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 8 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.002
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.250e-01

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.002
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.677e-13

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0008.vtu


## Step 10/41: t = 4.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 9 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00225
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.112e-01

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00225
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.439e-12

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0009.vtu


## Step 11/41: t = 5.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 10 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.0025
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.004e-01

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.0025
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.723e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0010.vtu


## Step 12/41: t = 5.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 11 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00275
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 9.183e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00275
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 9.336e-11

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0011.vtu


## Step 13/41: t = 6.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 12 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.003
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 8.539e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.003
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.670e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0012.vtu


## Step 14/41: t = 6.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 13 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00325
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 8.098e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00325
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 9.462e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0013.vtu


## Step 15/41: t = 7.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 14 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.0035
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 7.856e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.0035
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.768e-09

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0014.vtu


## Step 16/41: t = 7.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 15 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00375
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 7.789e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00375
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 8.190e-09

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0015.vtu


## Step 17/41: t = 8.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 16 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.004
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 7.823e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.004
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.539e-14

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0016.vtu


## Step 18/41: t = 8.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 17 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00425
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 7.807e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00425
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.007e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0017.vtu


## Step 19/41: t = 9.00e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 18 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.0045
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 7.539e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.0045
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.493e-09

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0018.vtu


## Step 20/41: t = 9.50e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 19 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00475
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 7.149e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00475
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.969e-09

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0019.vtu


## Step 21/41: t = 1.00e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 20 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.005
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 6.764e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.005
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.278e-09

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0020.vtu


## Step 22/41: t = 1.05e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 21 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00525
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 6.395e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00525
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.791e-09

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0021.vtu


## Step 23/41: t = 1.10e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 22 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.0055
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 6.056e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.0055
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.646e-09

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0022.vtu


## Step 24/41: t = 1.15e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 23 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00575
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.744e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00575
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.415e-09

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0023.vtu


## Step 25/41: t = 1.20e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 24 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.006
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.458e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.006
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.281e-09

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0024.vtu


## Step 26/41: t = 1.25e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 25 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00625
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.197e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00625
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.153e-09

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0025.vtu


## Step 27/41: t = 1.30e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 26 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.0065
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 4.958e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.0065
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 1.058e-09

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0026.vtu


## Step 28/41: t = 1.35e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 27 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00675
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 4.738e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00675
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 9.752e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0027.vtu


## Step 29/41: t = 1.40e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 28 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.007
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 4.535e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.007
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 9.067e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0028.vtu


## Step 30/41: t = 1.45e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 29 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00725
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 4.348e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00725
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 8.473e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0029.vtu


## Step 31/41: t = 1.50e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 30 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.0075
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 4.175e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.0075
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 7.960e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0030.vtu


## Step 32/41: t = 1.55e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 31 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00775
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 4.014e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00775
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 7.509e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0031.vtu


## Step 33/41: t = 1.60e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 32 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.008
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.865e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.008
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 7.109e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0032.vtu


## Step 34/41: t = 1.65e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 33 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00825
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.727e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00825
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 6.753e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0033.vtu


## Step 35/41: t = 1.70e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 34 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.0085
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.597e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.0085
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 6.432e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0034.vtu


## Step 36/41: t = 1.75e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 35 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00875
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.476e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00875
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 6.143e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0035.vtu


## Step 37/41: t = 1.80e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 36 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.009
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.362e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.009
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.880e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0036.vtu


## Step 38/41: t = 1.85e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 37 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00925
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.256e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00925
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.639e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0037.vtu


## Step 39/41: t = 1.90e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 38 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.0095
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.156e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.0095
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.419e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0038.vtu


## Step 40/41: t = 1.95e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 39 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00975
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 3.061e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.00975
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.217e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0039.vtu


## Step 41/41: t = 2.00e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 40 | dt = 5.00e-02 s
Coupling = staggered
  → Max iterations              : 50
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-04
  → Staggering tolerance |ΔD|   : 1.0e-04
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-07
  → Relative tolerance dmg      : 1.0e-06
  [Solver] Using quadrature degree 2 for integration measures.


#### Iteration 1/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.01
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 2.972e-02

Convergence check


#### Iteration 2/50


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 10 → [0.0, 0.0, 0.0]
  **[INFO]** Updating Displacement Dirichlet on region 11 → 0.01
  Building weak form, volume integrals (dx) for grain, tag = 1
  Non-linear solver (SNES Newton)
  ||Δu||/||u|| = 5.030e-10

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
[PlasticityModel] Updating plastic history...
Computing symbolic result fields (strain, stress, ...)
Exporting results to VTU file...
  → Projecting result fields for all materials...
Interpolation failed (likely Quadrature mismatch): Mismatch of tabulation points and element points.. Falling back to projection.
VTU file exported to: output/fields_0040.vtu

Simulation completed in 8.07 s
Total time steps solved: 41
