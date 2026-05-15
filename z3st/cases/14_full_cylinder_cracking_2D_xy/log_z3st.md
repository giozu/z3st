Info    : Reading 'mesh.msh'...
Info    : 15589 nodes
Info    : 31186 elements
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
  → Time steps          : 100
  → Regime              : 2d
  → Models active       :
      thermal    → ON
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
  → Temperature  : 0.8
  → Displacement : 0.6
  → Damage       : 0.5
  Adaptive relaxation enabled
  → relax_growth  : 1.1
  → relax_shrink : 0.95
  → relax_min  : 0.05
  → relax_max : 0.95


[ThermalModel] initializer
[ThermalModel] options loaded from input.yaml:
  analysis            : transient
  solver              : linear
  linear_solver       : iterative_hypre
  rtol                : 1e-06
  stag_tol            : 1e-06
  convergence         : rel_norm
[MechanicalModel] initializer
[MechanicalModel] options loaded from input.yaml:
  solver              : linear
  linear_solver       : direct_mumps
  rtol                : 1e-06
  stag_tol            : 1e-06
  convergence         : rel_norm
DamageModel initializer
Options loaded from input.yaml:
  type                : AT1
  solver              : linear
  linear_solver       : direct_mumps
  rtol                : 1e-06
  stag_tol            : 0.001
  convergence         : rel_norm
  lc                  : 5e-05
  hybrid_constraint   : True
[spine.load_materials]
Material loaded: uo2
  → k defined as constant: 5.0
  → Gc not defined for uo2
  - Material 'uo2': Gc (AT1) from sigma_c = 2.00e+09 Pa
  → constitutive model: lame
  E               → 358000000000.0 (float)
  G               → 145528455284.55286 (float)
  Gc              → 1489.75791433892 (float)
  T_initial       → 1023.15 (float)
  T_ref           → 298.15 (float)
  alpha           → 1e-05 (float)
  bulk_modulus    → 220987654320.98764 (float)
  constitutive_mode → lame (str)
  cp              → 280.0 (float)
  k               → 5.0 (float)
  lmbda           → 123968684131.28575 (float)
  name            → UO2 (str)
  nu              → 0.23 (float)
  rho             → 10970.0 (float)
  sigma_c         → 2000000000.0 (float)
[spine.initialize_fields]
[UPDATING q_third]

Initializing the temperature field...
  → Setting initial temperature for material: 'uo2'
    Set 15589 DOFs to 1023.15 K
  Initial T: min=1023.15 K, max=1023.15 K, mean=1023.15 K

Initializing the displacement field...
  Initial u: min=0.00e+00 m, max=0.00e+00 m, mean=0.00e+00 m

Initializing the damage field...



***


### spine - set_boundary_conditions


***



Loading boundary conditions from 'boundary_conditions.yaml'
  **[INFO]** Dirichlet thermal BC on 'uo2' → 263.15 K at region 'contact_wall'
  **[INFO]** Clamp_y mechanical BC on 'uo2' → 0.0 (first step) at region 'symmetry'
  **[INFO]** Clamp_x mechanical BC on 'uo2' → 0.0 (first step) at region 'pin'

Setting damage boundary conditions...
  **[INFO]** Dirichlet damage BC on 'uo2' → D = 1.0 at region 'crack_seed'


## Step 01/100: t = 1.00e-04 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 0 | dt = 1.00e-04 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=1011.42 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 8.650e-02
  [adaptive] relax_T=0.80

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.000e+00
  [adaptive] relax_u=0.60

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.000e+00
  [adaptive] relax_D=0.50
  |ΔD|_∞ = 5.000e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=1011.42 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.163e-03
  [adaptive] relax_T=0.88

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.000e-01
  [adaptive] relax_u=0.66

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.000e-01
  [adaptive] relax_D=0.55
  |ΔD|_∞ = 2.500e-01

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=1011.42 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.759e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.760e-01
  [adaptive] relax_u=0.73

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.116e-01
  [adaptive] relax_D=0.61
  |ΔD|_∞ = 2.183e-01

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=1011.42 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.165e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.582e-02
  [adaptive] relax_u=0.80

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.126e-01
  [adaptive] relax_D=0.67
  |ΔD|_∞ = 1.921e-01

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=1011.42 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.082e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.984e-02
  [adaptive] relax_u=0.88

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.494e-01
  [adaptive] relax_D=0.73
  |ΔD|_∞ = 1.059e-01

Convergence check


#### Iteration 6/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=1011.42 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.541e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.395e-03
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.104e-01
  [adaptive] relax_D=0.81
  |ΔD|_∞ = 4.542e-02

Convergence check


#### Iteration 7/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=1011.42 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.706e-09
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.777e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.535e-02
  [adaptive] relax_D=0.89
  |ΔD|_∞ = 1.431e-02

Convergence check


#### Iteration 8/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=1011.42 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.853e-10
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.888e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.731e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.115e-03

Convergence check


#### Iteration 9/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=1011.42 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.926e-11
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.444e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 9.556e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.843e-04

Convergence check


#### Iteration 10/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=1011.42 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.632e-13
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.221e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.821e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.935e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 10 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3588e+02 J
  → Fracture energy : 5.6374e-01 J
  → Total energy    : 2.3644e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0000.vtu


## Step 02/100: t = 1.11e-03 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 1 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=998.06 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.502e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.441e-03
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.977e-01
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.288e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=998.06 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.751e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.453e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.099e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.462e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=998.06 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.376e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.094e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.606e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.977e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=998.06 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.878e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.304e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.083e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.403e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=998.06 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.439e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.570e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 6.872e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.498e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3665e+02 J
  → Fracture energy : 5.3213e-01 J
  → Total energy    : 2.3718e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0001.vtu


## Step 03/100: t = 2.12e-03 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 2 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=990.02 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.585e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.016e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.504e-01
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.897e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=990.02 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.293e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.127e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.384e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.797e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=990.02 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.463e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.898e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 9.859e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.380e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=990.02 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.231e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.619e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 6.395e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.629e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=990.02 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.616e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.895e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.973e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.338e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3869e+02 J
  → Fracture energy : 5.2344e-01 J
  → Total energy    : 2.3921e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0002.vtu


## Step 04/100: t = 3.13e-03 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 3 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=984.02 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.737e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.304e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.208e-01
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.472e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=984.02 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.684e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.046e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.122e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.453e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=984.02 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.342e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.301e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 8.207e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.104e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=984.02 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.171e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.550e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.405e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.451e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=984.02 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.085e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.227e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.365e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.709e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4404e+02 J
  → Fracture energy : 5.1812e-01 J
  → Total energy    : 2.4455e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0003.vtu


## Step 05/100: t = 4.14e-03 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 4 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=979.10 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.332e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.739e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.032e-01
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.687e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=979.10 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.662e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.049e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.005e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.692e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=979.10 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.331e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.488e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.586e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.290e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=979.10 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.666e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.992e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.131e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.967e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=979.10 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.328e-08
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.874e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.266e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.042e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5929e+02 J
  → Fracture energy : 5.2005e-01 J
  → Total energy    : 2.5981e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0004.vtu


## Step 06/100: t = 5.15e-03 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 5 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=974.89 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.092e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.761e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 8.944e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.657e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=974.89 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.462e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.312e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 9.151e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.694e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=974.89 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.731e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.959e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.005e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.316e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=974.89 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.365e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.643e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.770e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.302e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=974.89 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.827e-08
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.655e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.051e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.182e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7899e+02 J
  → Fracture energy : 5.1981e-01 J
  → Total energy    : 2.7951e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0005.vtu


## Step 07/100: t = 6.15e-03 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 6 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=971.16 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.317e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.932e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.527e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.464e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=971.16 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.658e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.769e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.920e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.528e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=971.16 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.329e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.584e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 6.192e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.281e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=971.16 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.165e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.398e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.286e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.250e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=971.16 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.823e-08
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.503e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.777e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.181e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9472e+02 J
  → Fracture energy : 5.2420e-01 J
  → Total energy    : 2.9525e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0006.vtu


## Step 08/100: t = 7.16e-03 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 7 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=967.80 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 8.158e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.427e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 6.819e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.526e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=967.80 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.079e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.385e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.249e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.578e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=967.80 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.039e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.309e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.775e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.244e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=967.80 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.020e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.216e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.055e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.737e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=967.80 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.099e-08
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.390e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.655e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.750e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.0778e+02 J
  → Fracture energy : 5.3331e-01 J
  → Total energy    : 3.0831e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0007.vtu


## Step 09/100: t = 8.17e-03 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 8 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=964.72 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.279e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.102e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.204e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.431e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=964.72 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.640e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.094e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.396e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.465e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=964.72 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.820e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.093e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.787e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.157e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=964.72 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.099e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.073e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.024e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.098e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.2012e+02 J
  → Fracture energy : 5.4536e-01 J
  → Total energy    : 3.2067e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0008.vtu


## Step 10/100: t = 9.18e-03 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 9 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=961.88 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.587e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.837e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.156e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.501e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=961.88 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.294e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.848e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.370e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.517e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=961.88 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.647e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.911e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.737e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.174e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=961.88 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 8.234e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.952e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.966e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.121e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.3185e+02 J
  → Fracture energy : 5.6072e-01 J
  → Total energy    : 3.3241e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0009.vtu


## Step 11/100: t = 1.02e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 10 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=959.24 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.027e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.649e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.559e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.553e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=959.24 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.014e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.638e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.604e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.544e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=959.24 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.507e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.753e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.848e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.186e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=959.24 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.534e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.847e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.008e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.157e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.4404e+02 J
  → Fracture energy : 5.7914e-01 J
  → Total energy    : 3.4462e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0010.vtu


## Step 12/100: t = 1.12e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 11 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=956.75 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.563e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.558e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 8.078e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.577e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=956.75 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.782e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.456e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.825e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.547e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=956.75 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.391e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.610e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.913e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.171e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=956.75 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.954e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.751e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.007e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.946e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.5888e+02 J
  → Fracture energy : 6.0368e-01 J
  → Total energy    : 3.5948e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0011.vtu


## Step 13/100: t = 1.22e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 12 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=954.41 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.172e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.518e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 8.410e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.727e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=954.41 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.586e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.298e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 8.357e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.667e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=954.41 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.293e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.481e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 6.373e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.255e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=954.41 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.465e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.664e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.331e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.488e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.7692e+02 J
  → Fracture energy : 6.3246e-01 J
  → Total energy    : 3.7755e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0012.vtu


## Step 14/100: t = 1.32e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 13 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=952.19 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.838e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.493e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 8.766e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.840e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=952.19 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.419e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.158e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 8.510e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.772e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=952.19 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.209e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.365e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 6.437e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.329e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=952.19 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.047e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.586e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.353e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.959e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.9673e+02 J
  → Fracture energy : 6.6539e-01 J
  → Total energy    : 3.9740e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0013.vtu


## Step 15/100: t = 1.42e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 14 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=950.08 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.548e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.461e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 9.019e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.831e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=950.08 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.274e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.041e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 8.698e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.734e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=950.08 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.137e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.267e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 6.511e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.301e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=950.08 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.685e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.519e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.369e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.752e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.1720e+02 J
  → Fracture energy : 6.9791e-01 J
  → Total energy    : 4.1789e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0014.vtu


## Step 16/100: t = 1.52e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 15 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=948.07 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.294e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.541e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 8.681e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.915e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=948.07 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.147e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.959e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 8.481e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.824e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=948.07 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.073e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.183e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 6.407e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.362e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=948.07 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.367e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.459e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.323e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.154e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.3965e+02 J
  → Fracture energy : 7.3530e-01 J
  → Total energy    : 4.4038e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0015.vtu


## Step 17/100: t = 1.62e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 16 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=946.14 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.069e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.630e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 8.424e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.891e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=946.14 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.035e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.901e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 8.263e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.813e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=946.14 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.017e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.115e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 6.255e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.356e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=946.14 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.087e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.408e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.223e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.117e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.6519e+02 J
  → Fracture energy : 7.7233e-01 J
  → Total energy    : 4.6596e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0016.vtu


## Step 18/100: t = 1.73e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 17 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=944.29 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.869e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.718e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 8.326e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.908e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=944.29 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.935e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.853e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.940e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.799e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=944.29 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.674e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.054e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.939e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.340e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=944.29 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.837e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.362e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.988e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.002e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.9355e+02 J
  → Fracture energy : 8.0977e-01 J
  → Total energy    : 4.9436e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0017.vtu


## Step 19/100: t = 1.83e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 18 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=942.52 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.690e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.710e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 8.177e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.920e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=942.52 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.845e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.792e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.992e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.817e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=942.52 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.225e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.996e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 6.021e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.352e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=942.52 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.612e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.320e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.055e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.072e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.2293e+02 J
  → Fracture energy : 8.4664e-01 J
  → Total energy    : 5.2378e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0018.vtu


## Step 20/100: t = 1.93e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 19 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=940.80 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.528e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.730e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.960e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.903e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=940.80 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.764e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.743e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.804e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.817e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=940.80 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 8.820e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.945e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.915e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.357e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=940.80 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.410e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.282e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.004e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.120e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.5333e+02 J
  → Fracture energy : 8.8528e-01 J
  → Total energy    : 5.5421e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0019.vtu


## Step 21/100: t = 2.03e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 20 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=939.15 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.381e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.700e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.878e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.861e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=939.15 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.690e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.688e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.808e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.736e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=939.15 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 8.452e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.896e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.921e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.302e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=939.15 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.226e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.247e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.004e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.771e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.8482e+02 J
  → Fracture energy : 9.2334e-01 J
  → Total energy    : 5.8574e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0020.vtu


## Step 22/100: t = 2.13e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 21 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=937.55 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.247e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.666e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.859e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.829e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=937.55 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.623e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.634e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.829e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.732e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=937.55 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 8.117e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.849e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.970e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.291e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=937.55 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.058e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.214e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.053e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.677e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.1690e+02 J
  → Fracture energy : 9.6118e-01 J
  → Total energy    : 6.1786e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0021.vtu


## Step 23/100: t = 2.23e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 22 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=936.00 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.124e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.640e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.749e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.748e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=936.00 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.562e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.590e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.541e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.665e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=936.00 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.809e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.809e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.687e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.248e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=936.00 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.905e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.185e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.839e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.419e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.5009e+02 J
  → Fracture energy : 9.9780e-01 J
  → Total energy    : 6.5108e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0022.vtu


## Step 24/100: t = 2.33e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 23 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=934.50 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.010e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.595e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.324e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.716e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=934.50 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.505e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.542e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.119e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.634e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=934.50 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.526e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.770e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.354e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.220e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=934.50 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.763e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.157e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.608e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.214e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.8441e+02 J
  → Fracture energy : 1.0337e+00 J
  → Total energy    : 6.8544e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0023.vtu


## Step 25/100: t = 2.43e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 24 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=933.04 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.906e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.519e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 6.760e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.671e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=933.04 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.453e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.485e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 6.581e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.610e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=933.04 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.265e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.729e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.957e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.209e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=933.04 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.632e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.130e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.346e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.155e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.1903e+02 J
  → Fracture energy : 1.0701e+00 J
  → Total energy    : 7.2010e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0024.vtu


## Step 26/100: t = 2.53e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 25 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=931.63 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.809e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.443e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 6.120e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.641e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=931.63 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.404e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.434e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 6.068e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.550e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=931.63 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.022e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.692e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.632e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.155e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=931.63 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.511e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.105e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.151e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.775e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.5372e+02 J
  → Fracture energy : 1.1063e+00 J
  → Total energy    : 7.5483e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0025.vtu


## Step 27/100: t = 2.63e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 26 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=930.25 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.719e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.341e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.897e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.560e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=930.25 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.359e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.376e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.704e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.505e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=930.25 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.797e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.654e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.298e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.132e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=930.25 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.399e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.081e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.904e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.646e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.8853e+02 J
  → Fracture energy : 1.1402e+00 J
  → Total energy    : 7.8967e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0026.vtu


## Step 28/100: t = 2.73e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 27 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=928.91 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.635e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.281e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.439e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.521e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=928.91 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.317e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.330e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.482e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.478e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=928.91 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.587e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.621e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.199e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.113e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=928.91 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.294e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.059e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.862e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.530e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.2331e+02 J
  → Fracture energy : 1.1743e+00 J
  → Total energy    : 8.2448e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0027.vtu


## Step 29/100: t = 2.84e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 28 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=927.60 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.556e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.217e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.390e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.489e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=927.60 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.278e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.287e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.232e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.426e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=927.60 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.391e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.591e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.938e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.068e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=927.60 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.195e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.038e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.660e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.207e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.5834e+02 J
  → Fracture energy : 1.2076e+00 J
  → Total energy    : 8.5955e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0028.vtu


## Step 30/100: t = 2.94e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 29 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=926.33 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.483e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.160e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.005e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.432e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=926.33 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.241e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.247e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.963e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.392e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=926.33 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.207e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.562e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.790e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.050e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=926.33 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.103e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.019e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.582e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.108e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.9389e+02 J
  → Fracture energy : 1.2414e+00 J
  → Total energy    : 8.9513e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0029.vtu


## Step 31/100: t = 3.04e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 30 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=925.09 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.414e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.051e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.817e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.412e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=925.09 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.207e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.192e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.771e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.353e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=925.09 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.034e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.529e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.631e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.015e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=925.09 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.017e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.993e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.468e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.848e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.2929e+02 J
  → Fracture energy : 1.2731e+00 J
  → Total energy    : 9.3056e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0030.vtu


## Step 32/100: t = 3.14e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 31 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=923.87 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.349e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.019e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.635e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.334e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=923.87 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.174e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.162e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.568e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.300e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=923.87 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.871e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.505e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.463e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.826e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=923.87 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.936e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.823e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.348e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.663e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.6545e+02 J
  → Fracture energy : 1.3045e+00 J
  → Total energy    : 9.6676e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0031.vtu


## Step 33/100: t = 3.24e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 32 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=922.68 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.287e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.993e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.339e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.336e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=922.68 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.144e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.135e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.331e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.284e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=922.68 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.718e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.483e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.313e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.642e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=922.68 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.859e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.666e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.259e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.516e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0028e+03 J
  → Fracture energy : 1.3356e+00 J
  → Total energy    : 1.0041e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0032.vtu


## Step 34/100: t = 3.34e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 33 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=921.52 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.229e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.892e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.234e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.241e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=921.52 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.115e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.085e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.202e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.218e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=921.52 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.573e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.454e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.197e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.228e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=921.52 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.787e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.492e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.174e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.270e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0398e+03 J
  → Fracture energy : 1.3651e+00 J
  → Total energy    : 1.0412e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0033.vtu


## Step 35/100: t = 3.44e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 34 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=920.38 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.174e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.838e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.029e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.251e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=920.38 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.087e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.048e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.919e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.227e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=920.38 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.436e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.428e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.956e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.308e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=920.38 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.718e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.327e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.000e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.323e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0767e+03 J
  → Fracture energy : 1.3950e+00 J
  → Total energy    : 1.0781e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0034.vtu


## Step 36/100: t = 3.54e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 35 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=919.27 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.122e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.705e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.815e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.202e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=919.27 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.061e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.990e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.823e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.179e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=919.27 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.306e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.399e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.917e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.928e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=919.27 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.653e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.158e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.987e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.060e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1125e+03 J
  → Fracture energy : 1.4225e+00 J
  → Total energy    : 1.1139e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0035.vtu


## Step 37/100: t = 3.64e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 36 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=918.18 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.073e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.649e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.742e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.179e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=918.18 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.036e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.955e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.707e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.143e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=918.18 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.182e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.375e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.818e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.617e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=918.18 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.591e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.008e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.916e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.834e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1480e+03 J
  → Fracture energy : 1.4502e+00 J
  → Total energy    : 1.1495e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0036.vtu


## Step 38/100: t = 3.74e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 37 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=917.11 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.026e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.611e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.501e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.153e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=917.11 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.013e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.926e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.523e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.137e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=917.11 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.065e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.354e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.703e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.640e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=917.11 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.532e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.870e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.847e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.877e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1835e+03 J
  → Fracture energy : 1.4781e+00 J
  → Total energy    : 1.1850e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0037.vtu


## Step 39/100: t = 3.84e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 38 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=916.06 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.981e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.510e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.356e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.122e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=916.06 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 9.906e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.881e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.357e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.098e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=916.06 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.953e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.330e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.565e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.315e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=916.06 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.477e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.728e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.750e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.647e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2182e+03 J
  → Fracture energy : 1.5049e+00 J
  → Total energy    : 1.2197e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0038.vtu


## Step 40/100: t = 3.95e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 39 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=915.03 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.939e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.443e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.285e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.103e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=915.03 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 9.693e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.846e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.241e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.067e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=915.03 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.846e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.309e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.461e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.037e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=915.03 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.423e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.599e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.672e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.447e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2522e+03 J
  → Fracture energy : 1.5309e+00 J
  → Total energy    : 1.2537e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0039.vtu


## Step 41/100: t = 4.05e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 40 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=914.02 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.898e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.439e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.167e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.042e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=914.02 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 9.489e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.830e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.165e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.030e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=914.02 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.745e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.294e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.416e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.834e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=914.02 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.372e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.487e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.647e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.337e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2867e+03 J
  → Fracture energy : 1.5570e+00 J
  → Total energy    : 1.2883e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0040.vtu


## Step 42/100: t = 4.15e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 41 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=913.02 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.859e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.427e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.109e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.061e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=913.02 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 9.295e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.811e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.054e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.038e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=913.02 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.647e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.278e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.314e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.858e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=913.02 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.324e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.378e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.572e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.337e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3218e+03 J
  → Fracture energy : 1.5826e+00 J
  → Total energy    : 1.3234e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0041.vtu


## Step 43/100: t = 4.25e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 42 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=912.05 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.822e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.398e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.943e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.045e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=912.05 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 9.109e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.790e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.960e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.037e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=912.05 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.554e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.262e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.274e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.893e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=912.05 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.277e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.272e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.556e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.370e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3574e+03 J
  → Fracture energy : 1.6076e+00 J
  → Total energy    : 1.3590e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0042.vtu


## Step 44/100: t = 4.35e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 43 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=911.09 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.786e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.331e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.957e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.015e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=911.09 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.930e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.756e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.986e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.003e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=911.09 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.465e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.243e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.294e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.623e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=911.09 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.233e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.155e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.570e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.184e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3927e+03 J
  → Fracture energy : 1.6313e+00 J
  → Total energy    : 1.3943e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0043.vtu


## Step 45/100: t = 4.45e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 44 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=910.15 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.752e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.267e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.920e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.008e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=910.15 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.760e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.723e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.879e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.880e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=910.15 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.380e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.223e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.185e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.480e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=910.15 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.190e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.036e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.485e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.075e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4274e+03 J
  → Fracture energy : 1.6548e+00 J
  → Total energy    : 1.4291e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0044.vtu


## Step 46/100: t = 4.55e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 45 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=909.22 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.719e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.190e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.704e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.471e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=909.22 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.596e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.688e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.776e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.321e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=909.22 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.298e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.204e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.145e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.121e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=909.22 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.149e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.925e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.472e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.862e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4612e+03 J
  → Fracture energy : 1.6777e+00 J
  → Total energy    : 1.4629e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0045.vtu


## Step 47/100: t = 4.65e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 46 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=908.31 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.688e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.125e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.745e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.169e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=908.31 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.438e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.657e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.770e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.100e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=908.31 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.219e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.187e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.125e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.938e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=908.31 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.110e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.820e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.453e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.733e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4942e+03 J
  → Fracture energy : 1.6994e+00 J
  → Total energy    : 1.4959e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0046.vtu


## Step 48/100: t = 4.75e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 47 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=907.41 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.657e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.124e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.702e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.148e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=907.41 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.287e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.643e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.681e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.010e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=907.41 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.143e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.173e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.052e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.845e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=907.41 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.072e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.728e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.403e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.669e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5276e+03 J
  → Fracture energy : 1.7216e+00 J
  → Total energy    : 1.5293e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0047.vtu


## Step 49/100: t = 4.85e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 48 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=906.53 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.628e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.079e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.532e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.870e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=906.53 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.141e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.618e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.581e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.849e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=906.53 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.070e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.158e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.994e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.755e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=906.53 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.035e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.633e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.370e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.609e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5607e+03 J
  → Fracture energy : 1.7432e+00 J
  → Total energy    : 1.5624e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0048.vtu


## Step 50/100: t = 4.95e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 49 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=905.66 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.600e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.012e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.505e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.660e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=905.66 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.000e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.590e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.505e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.594e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=905.66 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.000e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.143e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.915e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.547e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=905.66 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.000e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.541e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.307e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.463e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5930e+03 J
  → Fracture energy : 1.7638e+00 J
  → Total energy    : 1.5948e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0049.vtu


## Step 51/100: t = 5.06e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 50 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=904.80 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.573e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.981e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.404e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.906e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=904.80 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.865e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.570e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.407e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.752e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=904.80 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.933e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.129e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.844e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.636e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=904.80 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.966e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.457e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.261e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.508e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6251e+03 J
  → Fracture energy : 1.7842e+00 J
  → Total energy    : 1.6269e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0050.vtu


## Step 52/100: t = 5.16e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 51 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=903.96 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.547e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.943e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.336e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.628e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=903.96 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.735e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.549e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.370e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.402e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=903.96 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.867e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.116e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.822e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.345e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=903.96 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.934e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.372e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.248e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.300e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6568e+03 J
  → Fracture energy : 1.8036e+00 J
  → Total energy    : 1.6586e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0051.vtu


## Step 53/100: t = 5.26e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 52 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=903.13 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.522e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.932e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.313e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.320e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=903.13 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.609e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.533e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.317e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.253e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=903.13 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.804e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.104e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.778e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.272e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=903.13 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.902e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.290e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.216e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.265e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6887e+03 J
  → Fracture energy : 1.8230e+00 J
  → Total energy    : 1.6905e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0052.vtu


## Step 54/100: t = 5.36e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 53 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=902.31 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.497e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.892e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.311e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.232e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=902.31 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.487e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.512e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.314e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.127e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=902.31 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.744e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.091e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.768e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.163e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=902.31 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.872e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.210e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.206e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.181e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7202e+03 J
  → Fracture energy : 1.8418e+00 J
  → Total energy    : 1.7220e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0053.vtu


## Step 55/100: t = 5.46e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 54 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=901.50 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.474e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.880e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.256e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.455e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=901.50 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.370e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.498e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.233e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.445e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=901.50 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.685e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.080e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.697e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.444e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=901.50 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.842e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.133e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.153e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.393e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7518e+03 J
  → Fracture energy : 1.8608e+00 J
  → Total energy    : 1.7536e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0054.vtu


## Step 56/100: t = 5.56e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 55 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=900.70 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.451e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.844e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.095e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.304e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=900.70 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.256e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.479e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.135e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.274e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=900.70 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.628e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.068e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.643e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.309e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=900.70 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.814e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.060e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.125e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.299e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7831e+03 J
  → Fracture energy : 1.8794e+00 J
  → Total energy    : 1.7849e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0055.vtu


## Step 57/100: t = 5.66e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 56 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=899.91 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.429e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.798e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.048e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.971e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=899.91 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.146e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.458e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 2.053e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.923e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=899.91 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.573e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.056e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.568e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.037e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=899.91 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.787e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.989e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.070e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.113e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8139e+03 J
  → Fracture energy : 1.8974e+00 J
  → Total energy    : 1.8158e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0056.vtu


## Step 58/100: t = 5.76e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 57 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=899.14 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.408e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.759e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.948e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.689e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=899.14 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.040e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.439e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.973e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.607e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=899.14 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.520e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.045e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.514e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.786e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=899.14 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.760e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.918e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.035e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.938e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8442e+03 J
  → Fracture energy : 1.9150e+00 J
  → Total energy    : 1.8461e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0057.vtu


## Step 59/100: t = 5.86e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 58 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=898.37 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.387e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.723e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.909e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.251e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=898.37 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.937e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.420e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.910e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.223e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=898.37 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.468e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.033e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.459e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.542e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=898.37 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.734e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.848e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 9.951e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.793e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8741e+03 J
  → Fracture energy : 1.9320e+00 J
  → Total energy    : 1.8761e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0058.vtu


## Step 60/100: t = 5.96e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 59 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=897.61 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.367e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.709e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.807e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.344e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=897.61 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.837e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.406e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.795e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.365e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=897.61 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.419e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.023e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.366e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.636e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=897.61 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.709e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.781e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 9.300e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.852e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9040e+03 J
  → Fracture energy : 1.9491e+00 J
  → Total energy    : 1.9060e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0059.vtu


## Step 61/100: t = 6.06e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 60 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=896.87 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.348e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.680e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.639e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.184e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=896.87 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.740e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.390e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.657e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.182e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=896.87 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.370e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.013e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.277e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.491e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=896.87 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.685e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.718e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 8.757e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.750e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9336e+03 J
  → Fracture energy : 1.9661e+00 J
  → Total energy    : 1.9355e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0060.vtu


## Step 62/100: t = 6.17e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 61 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=896.13 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.329e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.639e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.638e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.837e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=896.13 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.646e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.372e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.653e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.824e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=896.13 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.323e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.003e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.267e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.215e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=896.13 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.662e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.656e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 8.659e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.562e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9626e+03 J
  → Fracture energy : 1.9824e+00 J
  → Total energy    : 1.9646e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0061.vtu


## Step 63/100: t = 6.27e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 62 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=895.40 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.311e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.649e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.649e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.046e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=895.40 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.555e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.364e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.642e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.982e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=895.40 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.277e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.938e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.251e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.315e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=895.40 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.639e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.592e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 8.512e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.620e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9922e+03 J
  → Fracture energy : 1.9987e+00 J
  → Total energy    : 1.9942e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0062.vtu


## Step 64/100: t = 6.37e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 63 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=894.68 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.293e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.639e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.573e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.027e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=894.68 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.466e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.352e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.593e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.065e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=894.68 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.233e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.849e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.220e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.400e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=894.68 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.617e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.532e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 8.328e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.682e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0218e+03 J
  → Fracture energy : 2.0150e+00 J
  → Total energy    : 2.0238e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0063.vtu


## Step 65/100: t = 6.47e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 64 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=893.97 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.276e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.608e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.516e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.031e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=893.97 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.380e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.337e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.509e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.053e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=893.97 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.190e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.756e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.147e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.386e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=893.97 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.595e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.475e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.795e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.671e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0510e+03 J
  → Fracture energy : 2.0308e+00 J
  → Total energy    : 2.0530e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0064.vtu


## Step 66/100: t = 6.57e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 65 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=893.27 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.259e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.567e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.422e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.922e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=893.27 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.297e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.320e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.451e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.927e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=893.27 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.148e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.667e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.119e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.286e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=893.27 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.574e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.422e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.665e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.601e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0795e+03 J
  → Fracture energy : 2.0460e+00 J
  → Total energy    : 2.0816e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0065.vtu


## Step 67/100: t = 6.67e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 66 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=892.57 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.243e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.534e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.417e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.771e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=892.57 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.215e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.305e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.447e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.755e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=892.57 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.108e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.578e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.114e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.150e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=892.57 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.554e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.368e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.629e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.507e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1076e+03 J
  → Fracture energy : 2.0609e+00 J
  → Total energy    : 2.1096e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0066.vtu


## Step 68/100: t = 6.77e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 67 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=891.89 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.227e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.530e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.457e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.018e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=891.89 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.136e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.296e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.451e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.956e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=891.89 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.068e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.500e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.105e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.284e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=891.89 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.534e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.314e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.516e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.590e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1357e+03 J
  → Fracture energy : 2.0760e+00 J
  → Total energy    : 2.1378e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0067.vtu


## Step 69/100: t = 6.87e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 68 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=891.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.212e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.537e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.423e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.155e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=891.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.059e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.288e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.430e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.039e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=891.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.029e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.425e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.094e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.328e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=891.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.515e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.261e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.466e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.611e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1642e+03 J
  → Fracture energy : 2.0914e+00 J
  → Total energy    : 2.1663e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0068.vtu


## Step 70/100: t = 6.97e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 69 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=890.54 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.197e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.525e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.381e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.905e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=890.54 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.984e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.278e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.373e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.753e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=890.54 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.992e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.346e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.045e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.129e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=890.54 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.496e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.209e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.105e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.477e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1928e+03 J
  → Fracture energy : 2.1068e+00 J
  → Total energy    : 2.1949e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0069.vtu


## Step 71/100: t = 7.07e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 70 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=889.87 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.182e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.499e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.335e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.510e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=889.87 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.911e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.264e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.366e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.494e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=889.87 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.955e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.261e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.057e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.934e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=889.87 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.478e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.156e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.258e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.346e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2210e+03 J
  → Fracture energy : 2.1215e+00 J
  → Total energy    : 2.2231e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0070.vtu


## Step 72/100: t = 7.17e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 71 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=889.22 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.168e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.483e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.369e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.287e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=889.22 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.840e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.253e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.364e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.268e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=889.22 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.920e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.187e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.037e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.762e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=889.22 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.460e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.109e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.037e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.230e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2491e+03 J
  → Fracture energy : 2.1359e+00 J
  → Total energy    : 2.2512e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0071.vtu


## Step 73/100: t = 7.28e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 72 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=888.57 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.154e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.475e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.337e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.336e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=888.57 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.770e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.244e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.377e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.304e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=888.57 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.885e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.116e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.065e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.785e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=888.57 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.443e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.061e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.299e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.244e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2772e+03 J
  → Fracture energy : 2.1502e+00 J
  → Total energy    : 2.2793e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0072.vtu


## Step 74/100: t = 7.38e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 73 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=887.93 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.141e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.454e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.372e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.075e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=887.93 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.703e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.231e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.373e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.038e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=887.93 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.851e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.040e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.046e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.583e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=887.93 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.426e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.013e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.108e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.107e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3050e+03 J
  → Fracture energy : 2.1641e+00 J
  → Total energy    : 2.3072e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0073.vtu


## Step 75/100: t = 7.48e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 74 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=887.29 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.127e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.439e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.320e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.962e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=887.29 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.637e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.221e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.321e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.908e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=887.29 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.818e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.968e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.005e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.525e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=887.29 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.409e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.967e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 6.831e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.091e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3326e+03 J
  → Fracture energy : 2.1780e+00 J
  → Total energy    : 2.3348e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0074.vtu


## Step 76/100: t = 7.58e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 75 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=886.66 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.114e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.405e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.219e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.748e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=886.66 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.572e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.208e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.219e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.740e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=886.66 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.786e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.897e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 9.285e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.372e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=886.66 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.393e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.925e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 6.313e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.989e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3596e+03 J
  → Fracture energy : 2.1916e+00 J
  → Total energy    : 2.3618e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0075.vtu


## Step 77/100: t = 7.68e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 76 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=886.04 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.102e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.361e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.098e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.469e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=886.04 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.509e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.194e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.099e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.462e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=886.04 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.755e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.830e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 8.372e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.139e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=886.04 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.377e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.887e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.696e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.796e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3857e+03 J
  → Fracture energy : 2.2045e+00 J
  → Total energy    : 2.3879e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0076.vtu


## Step 78/100: t = 7.78e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 77 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=885.42 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.090e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.315e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.036e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.083e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=885.42 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.448e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.181e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.069e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.080e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=885.42 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.724e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.771e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 8.267e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.854e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=885.42 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.362e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.854e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.677e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.607e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4106e+03 J
  → Fracture energy : 2.2167e+00 J
  → Total energy    : 2.4128e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0077.vtu


## Step 79/100: t = 7.88e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 78 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=884.81 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.078e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.296e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.069e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.014e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=884.81 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.388e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.171e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.089e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.006e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=884.81 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.694e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.710e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 8.365e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.796e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=884.81 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.347e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.816e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.713e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.568e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4351e+03 J
  → Fracture energy : 2.2285e+00 J
  → Total energy    : 2.4373e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0078.vtu


## Step 80/100: t = 7.98e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 79 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=884.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.066e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.297e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.097e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.016e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=884.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.330e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.163e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.099e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.071e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=884.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.665e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.646e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 8.373e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.879e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=884.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.332e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.772e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.695e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.644e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4598e+03 J
  → Fracture energy : 2.2406e+00 J
  → Total energy    : 2.4620e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0079.vtu


## Step 81/100: t = 8.08e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 80 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=883.61 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.054e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.297e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.078e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.064e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=883.61 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.272e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.156e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.104e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.118e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=883.61 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.636e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.584e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 8.483e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.915e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=883.61 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.318e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.730e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.792e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.668e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4847e+03 J
  → Fracture energy : 2.2529e+00 J
  → Total energy    : 2.4870e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0080.vtu


## Step 82/100: t = 8.18e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 81 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=883.02 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.043e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.289e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.099e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.123e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=883.02 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.216e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.148e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.100e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.172e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=883.02 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.608e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.522e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 8.371e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.955e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=883.02 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.304e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.689e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.683e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.695e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5096e+03 J
  → Fracture energy : 2.2652e+00 J
  → Total energy    : 2.5119e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0081.vtu


## Step 83/100: t = 8.28e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 82 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=882.43 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.032e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.278e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.060e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.178e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=882.43 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.162e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.139e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.061e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.220e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=882.43 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.581e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.461e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 8.073e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.990e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=882.43 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.290e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.650e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.478e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.718e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5345e+03 J
  → Fracture energy : 2.2774e+00 J
  → Total energy    : 2.5367e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0082.vtu


## Step 84/100: t = 8.39e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 83 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=881.85 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.022e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.253e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 9.858e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.026e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=881.85 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.108e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.128e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 9.875e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.067e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=881.85 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.554e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.401e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.513e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.875e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=881.85 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.277e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.614e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.101e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.642e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5588e+03 J
  → Fracture energy : 2.2892e+00 J
  → Total energy    : 2.5611e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0083.vtu


## Step 85/100: t = 8.49e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 84 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=881.27 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.011e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.222e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 9.448e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.728e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=881.27 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 5.056e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.118e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 9.750e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.772e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=881.27 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.528e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.346e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.578e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.654e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=881.27 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.264e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.582e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.219e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.495e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5823e+03 J
  → Fracture energy : 2.3004e+00 J
  → Total energy    : 2.5847e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0084.vtu


## Step 86/100: t = 8.59e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 85 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=880.70 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.001e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.202e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 9.901e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.742e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=880.70 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 5.005e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.107e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 9.980e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.776e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=880.70 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.502e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.284e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.635e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.654e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=880.70 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.251e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.544e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.208e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.493e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6055e+03 J
  → Fracture energy : 2.3112e+00 J
  → Total energy    : 2.6078e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0085.vtu


## Step 87/100: t = 8.69e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 86 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=880.14 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 9.909e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.197e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.006e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.956e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=880.14 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.955e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.098e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.013e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.966e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=880.14 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.477e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.213e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.739e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.789e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=880.14 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.239e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.498e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.274e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.580e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6289e+03 J
  → Fracture energy : 2.3224e+00 J
  → Total energy    : 2.6312e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0086.vtu


## Step 88/100: t = 8.79e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 87 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=879.58 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 9.811e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.192e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.003e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.061e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=879.58 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.906e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.089e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.007e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.046e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=879.58 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.453e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.150e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.684e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.840e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=879.58 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.226e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.457e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.233e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.610e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6523e+03 J
  → Fracture energy : 2.3340e+00 J
  → Total energy    : 2.6547e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0087.vtu


## Step 89/100: t = 8.89e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 88 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=879.02 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 9.715e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.193e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 9.942e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.034e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=879.02 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.858e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.083e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.010e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.997e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=879.02 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.429e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.092e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.753e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.795e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=879.02 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.214e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.419e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.300e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.575e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.6760e+03 J
  → Fracture energy : 2.3459e+00 J
  → Total energy    : 2.6783e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0088.vtu


## Step 90/100: t = 8.99e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 89 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=878.47 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 9.621e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.197e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.000e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.873e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=878.47 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.811e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.076e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 1.002e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.816e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=878.47 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.405e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.035e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.652e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.652e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=878.47 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.203e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.379e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.212e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.495e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7000e+03 J
  → Fracture energy : 2.3579e+00 J
  → Total energy    : 2.7023e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0089.vtu


## Step 91/100: t = 9.09e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 90 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=877.93 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 9.529e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.191e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 9.647e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.703e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=877.93 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.764e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.069e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 9.788e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.753e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=877.93 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.382e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.979e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.497e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.644e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=877.93 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.191e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.343e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.111e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.490e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7240e+03 J
  → Fracture energy : 2.3695e+00 J
  → Total energy    : 2.7264e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0090.vtu


## Step 92/100: t = 9.19e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 91 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=877.39 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 9.439e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.176e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 9.757e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.465e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=877.39 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.719e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.061e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 9.903e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.511e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=877.39 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.360e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.932e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.574e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.463e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=877.39 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.180e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.313e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.157e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.369e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7478e+03 J
  → Fracture energy : 2.3807e+00 J
  → Total energy    : 2.7502e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0091.vtu


## Step 93/100: t = 9.29e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 92 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=876.85 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 9.350e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.173e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 9.788e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.597e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=876.85 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.675e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.056e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 9.798e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.594e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=876.85 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.337e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.885e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.449e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.485e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=876.85 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.169e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.281e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.053e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.357e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7716e+03 J
  → Fracture energy : 2.3917e+00 J
  → Total energy    : 2.7740e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0092.vtu


## Step 94/100: t = 9.39e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 93 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=876.32 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 9.263e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.160e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 9.296e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.536e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=876.32 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.632e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.049e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 9.466e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.531e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=876.32 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.316e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.844e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.271e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.437e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=876.32 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.158e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.254e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.963e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.325e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7952e+03 J
  → Fracture energy : 2.4026e+00 J
  → Total energy    : 2.7976e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0093.vtu


## Step 95/100: t = 9.50e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 94 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=875.80 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 9.178e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.138e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 9.440e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.285e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=875.80 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.589e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.042e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 9.579e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.283e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=875.80 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.295e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.804e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.318e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.269e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=875.80 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.147e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.230e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.976e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.227e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8181e+03 J
  → Fracture energy : 2.4129e+00 J
  → Total energy    : 2.8205e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0094.vtu


## Step 96/100: t = 9.60e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 95 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=875.27 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 9.095e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.132e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 9.545e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.416e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=875.27 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.547e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.035e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 9.579e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.465e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=875.27 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.274e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.754e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.287e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.413e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=875.27 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.137e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.197e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.944e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.324e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8410e+03 J
  → Fracture energy : 2.4232e+00 J
  → Total energy    : 2.8435e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0095.vtu


## Step 97/100: t = 9.70e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 96 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=874.76 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 9.013e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.121e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 9.186e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.440e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=874.76 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.506e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.028e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 9.221e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.488e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=874.76 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.253e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.709e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.017e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.431e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=874.76 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.127e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.169e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.763e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.337e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8638e+03 J
  → Fracture energy : 2.4336e+00 J
  → Total energy    : 2.8662e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0096.vtu


## Step 98/100: t = 9.80e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 97 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=874.24 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 8.932e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.106e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 8.626e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.328e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=874.24 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.466e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.021e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 8.666e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.377e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=874.24 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.233e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.668e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 6.601e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.349e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=874.24 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.117e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.143e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.484e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.283e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.8862e+03 J
  → Fracture energy : 2.4439e+00 J
  → Total energy    : 2.8886e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0097.vtu


## Step 99/100: t = 9.90e-02 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 98 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=873.73 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 8.853e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.086e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.949e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.129e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=873.73 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.427e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.015e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.998e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.180e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=873.73 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.213e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.632e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 6.100e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.203e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=873.73 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.107e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.121e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 4.149e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.186e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9080e+03 J
  → Fracture energy : 2.4537e+00 J
  → Total energy    : 2.9104e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0098.vtu


## Step 100/100: t = 1.00e-01 s | LHR = 0.00e+00 W/m

[UPDATING q_third]



***


### spine - solve


***



Current step = 99 | dt = 1.01e-03 s
Coupling = staggered
  → Max iterations              : 200
  → Staggering tolerance |ΔT|   : 1.0e-06
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-03
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=873.23 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 8.776e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.064e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.210e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.870e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=873.23 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.388e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.009e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 7.269e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.926e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=873.23 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.194e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.602e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 5.555e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.014e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=873.23 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.097e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.102e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=1.49e+03, sigma_c=2.00e+09
  ||ΔD||/||D|| = 3.786e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.061e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9291e+03 J
  → Fracture energy : 2.4631e+00 J
  → Total energy    : 2.9315e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0099.vtu

Simulation completed in 288.25 s
Total time steps solved: 100
