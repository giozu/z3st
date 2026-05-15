Info    : Reading 'mesh.msh'...
Info    : 17965 nodes
Info    : 35948 elements
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
  - Material 'uo2': Gc (AT1) from sigma_c = 1.00e+09 Pa
  → constitutive model: lame
  E               → 358000000000.0 (float)
  G               → 145528455284.55286 (float)
  Gc              → 372.43947858473 (float)
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
  sigma_c         → 1000000000.0 (float)
[spine.initialize_fields]
[UPDATING q_third]

Initializing the temperature field...
  → Setting initial temperature for material: 'uo2'
    Set 17965 DOFs to 1023.15 K
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
  T_new: min=263.15 K, max=1023.15 K, mean=996.32 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.170e-01
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.000e+00
  [adaptive] relax_D=0.50
  |ΔD|_∞ = 5.000e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=996.32 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.912e-03
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.459e-01
  [adaptive] relax_D=0.55
  |ΔD|_∞ = 3.348e-01

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=996.32 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.741e-03
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.658e+00
  [adaptive] relax_D=0.52
  |ΔD|_∞ = 4.326e-01

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=996.32 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.255e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.581e-02
  [adaptive] relax_u=0.80

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.706e+00
  [adaptive] relax_D=0.50
  |ΔD|_∞ = 2.138e-01

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=996.32 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.128e-05
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.348e+00
  [adaptive] relax_D=0.47
  |ΔD|_∞ = 1.087e-01

Convergence check


#### Iteration 6/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=996.32 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.638e-07
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.756e-01
  [adaptive] relax_D=0.52
  |ΔD|_∞ = 5.763e-02

Convergence check


#### Iteration 7/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=996.32 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.819e-08
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.776e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.702e-01
  [adaptive] relax_D=0.57
  |ΔD|_∞ = 3.430e-02

Convergence check


#### Iteration 8/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=996.32 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.409e-09
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.500e-01
  [adaptive] relax_D=0.63
  |ΔD|_∞ = 1.820e-02

Convergence check


#### Iteration 9/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=996.32 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.047e-11
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.182e-01
  [adaptive] relax_D=0.69
  |ΔD|_∞ = 8.600e-03

Convergence check


#### Iteration 10/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=996.32 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.524e-12
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.220e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.840e-02
  [adaptive] relax_D=0.76
  |ΔD|_∞ = 3.523e-03

Convergence check


#### Iteration 11/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=996.32 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.762e-13
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.610e-09
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.648e-02
  [adaptive] relax_D=0.84
  |ΔD|_∞ = 1.200e-03

Convergence check


#### Iteration 12/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=996.32 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 8.809e-15
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.805e-10
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.361e-03
  [adaptive] relax_D=0.92
  |ΔD|_∞ = 3.174e-04

Convergence check


#### Iteration 13/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=996.32 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.407e-16
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.025e-12
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.897e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.748e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 13 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3571e+02 J
  → Fracture energy : 3.1642e-01 J
  → Total energy    : 2.3603e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=959.68 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.208e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.397e-03
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.939e-01
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.875e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=959.68 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.604e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.424e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.469e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.441e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=959.68 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.302e-04
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.831e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.173e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=959.68 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.151e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.279e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.738e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.543e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=959.68 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.755e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.573e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.077e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.933e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.4070e+02 J
  → Fracture energy : 7.6716e-01 J
  → Total energy    : 2.4146e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=939.43 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.289e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.904e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.543e-01
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.294e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=939.43 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.145e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.447e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.891e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.985e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=939.43 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.072e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.227e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.137e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.887e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=939.43 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.362e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.667e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.874e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.537e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=939.43 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.681e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.369e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.382e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.137e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5678e+02 J
  → Fracture energy : 1.0278e+00 J
  → Total energy    : 2.5781e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=925.05 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.848e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.655e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.512e-01
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.744e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=925.05 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.424e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.362e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.101e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.550e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=925.05 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.120e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.345e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.515e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.444e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=925.05 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.560e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.018e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.068e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.162e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=925.05 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.780e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.934e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.692e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.536e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7706e+02 J
  → Fracture energy : 1.2138e+00 J
  → Total energy    : 2.7827e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=913.72 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.159e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.274e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.072e-01
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.479e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=913.72 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.079e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.881e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.763e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.374e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=913.72 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.396e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.886e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.659e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.240e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=913.72 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.698e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.673e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.010e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.749e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=913.72 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.349e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.703e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.023e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.057e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.9814e+02 J
  → Fracture energy : 1.3646e+00 J
  → Total energy    : 2.9950e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=904.30 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.751e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.084e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.343e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.368e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=904.30 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 8.754e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.537e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.089e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.282e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=904.30 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.377e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.570e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.766e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.125e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=904.30 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.188e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.441e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.492e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.688e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=904.30 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.094e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.549e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.679e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.195e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1888e+02 J
  → Fracture energy : 1.4944e+00 J
  → Total energy    : 3.2038e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=896.19 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.479e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.048e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.846e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.341e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=896.19 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.397e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.309e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.081e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.220e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=896.19 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.698e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.350e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.267e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.052e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=896.19 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.849e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.278e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.200e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.036e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=896.19 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.246e-08
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.442e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.480e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.679e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.4002e+02 J
  → Fracture energy : 1.6090e+00 J
  → Total energy    : 3.4163e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=889.04 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.285e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.123e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.828e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.383e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=889.04 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.424e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.166e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.438e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.164e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=889.04 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.212e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.195e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.970e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.952e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=889.04 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.606e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.160e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.019e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.516e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=889.04 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 8.030e-08
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.363e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.350e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.274e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.6279e+02 J
  → Fracture energy : 1.7113e+00 J
  → Total energy    : 3.6450e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=882.65 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.138e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.344e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.145e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.491e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=882.65 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.689e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.097e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.045e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.235e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=882.65 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.845e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.088e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.785e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.382e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=882.65 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.422e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.074e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.896e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.062e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=882.65 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.111e-08
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.303e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.258e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.940e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.8786e+02 J
  → Fracture energy : 1.8041e+00 J
  → Total energy    : 3.8966e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=876.85 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.023e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.772e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.688e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.617e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=876.85 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.113e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.104e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.799e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.360e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=876.85 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.557e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.019e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.665e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.232e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=876.85 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.278e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.008e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.806e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.629e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=876.85 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.391e-08
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.257e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.187e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.622e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.1622e+02 J
  → Fracture energy : 1.8904e+00 J
  → Total energy    : 4.1811e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=871.54 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.297e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.013e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.341e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.693e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=871.54 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.649e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.069e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.630e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.463e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=871.54 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.324e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.942e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.571e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.990e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=871.54 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.162e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.944e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.730e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.245e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=871.54 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.811e-08
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.213e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.126e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.302e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.5130e+02 J
  → Fracture energy : 1.9711e+00 J
  → Total energy    : 4.5327e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=866.64 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 8.531e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.875e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.134e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.711e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=866.64 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.266e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.924e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.528e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.494e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=866.64 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.133e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.832e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.500e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.044e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=866.64 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.066e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.871e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.668e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.628e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=866.64 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.332e-08
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.167e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.076e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.993e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.8693e+02 J
  → Fracture energy : 2.0484e+00 J
  → Total energy    : 4.8898e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=862.08 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.887e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.970e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.957e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.683e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=862.08 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.944e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.826e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.403e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.496e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=862.08 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.972e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.731e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.404e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.052e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=862.08 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.859e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.799e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.591e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.762e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.2435e+02 J
  → Fracture energy : 2.1239e+00 J
  → Total energy    : 5.2648e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=857.83 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.338e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.846e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.721e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.613e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=857.83 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.669e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.639e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.208e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.448e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=857.83 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.835e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.593e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.263e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.032e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=857.83 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.173e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.711e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.490e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.649e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.6333e+02 J
  → Fracture energy : 2.1953e+00 J
  → Total energy    : 5.6553e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=853.84 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.864e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.812e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.462e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.722e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=853.84 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.432e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.485e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.995e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.524e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=853.84 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.716e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.467e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.111e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.079e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=853.84 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 8.580e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.628e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.386e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.924e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.0470e+02 J
  → Fracture energy : 2.2614e+00 J
  → Total energy    : 6.0696e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=850.07 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.450e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.381e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.235e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.795e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=850.07 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.225e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.229e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.810e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.592e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=850.07 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.612e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.318e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.979e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.137e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=850.07 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 8.062e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.541e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.293e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.397e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.4495e+02 J
  → Fracture energy : 2.3222e+00 J
  → Total energy    : 6.4727e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=846.51 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.085e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.019e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.025e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.813e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=846.51 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.042e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.019e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.635e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.633e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=846.51 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.521e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.193e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.853e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.174e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=846.51 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.606e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.467e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.206e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.651e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 6.8234e+02 J
  → Fracture energy : 2.3787e+00 J
  → Total energy    : 6.8472e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=843.14 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.761e-03
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.798e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.765e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=843.14 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.880e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.844e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.449e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.600e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=843.14 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.440e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.084e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.724e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.162e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=843.14 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.201e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.400e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.120e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.653e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.1620e+02 J
  → Fracture energy : 2.4309e+00 J
  → Total energy    : 7.1863e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=839.92 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.471e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.514e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.571e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.646e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=839.92 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.735e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.708e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.257e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.506e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=839.92 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.368e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.995e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.589e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.104e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=839.92 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.838e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.342e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.031e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.339e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.4653e+02 J
  → Fracture energy : 2.4789e+00 J
  → Total energy    : 7.4901e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=836.86 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.209e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.347e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.338e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.567e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=836.86 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.605e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.603e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.060e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.394e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=836.86 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.302e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.922e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.451e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.025e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=836.86 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.512e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.294e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.402e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.874e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.7553e+02 J
  → Fracture energy : 2.5229e+00 J
  → Total energy    : 7.7805e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=833.93 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.973e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.107e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.121e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.514e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=833.93 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.487e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.490e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.871e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.357e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=833.93 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.243e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.853e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.318e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.611e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=833.93 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.216e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.251e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.545e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.267e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.0261e+02 J
  → Fracture energy : 2.5632e+00 J
  → Total energy    : 8.0518e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=831.12 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.758e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.893e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.914e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.439e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=831.12 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.379e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.391e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.691e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.300e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=831.12 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.189e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.792e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.192e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.248e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=831.12 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.947e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.211e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.743e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.929e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.2806e+02 J
  → Fracture energy : 2.6004e+00 J
  → Total energy    : 8.3066e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=828.43 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.561e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.710e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.727e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.343e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=828.43 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.281e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.305e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.530e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.223e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=828.43 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.140e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.738e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.080e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.729e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=828.43 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.702e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.176e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.018e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.656e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.5251e+02 J
  → Fracture energy : 2.6349e+00 J
  → Total energy    : 8.5514e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=825.84 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.381e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.535e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.562e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.269e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=825.84 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.190e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.228e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.386e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.163e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=825.84 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.095e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.690e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.807e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.315e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=825.84 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.476e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.145e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.387e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.380e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.7584e+02 J
  → Fracture energy : 2.6667e+00 J
  → Total energy    : 8.7851e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=823.35 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.214e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.391e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.410e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.198e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=823.35 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.107e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.167e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.257e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.103e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=823.35 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.054e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.650e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.934e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.010e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=823.35 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.268e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.119e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.838e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.246e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.9822e+02 J
  → Fracture energy : 2.6965e+00 J
  → Total energy    : 9.0091e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=820.94 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.061e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.284e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.279e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.110e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=820.94 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.030e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.120e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.147e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.031e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=820.94 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.015e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.618e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.196e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.478e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=820.94 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.076e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.096e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.374e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.880e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.2007e+02 J
  → Fracture energy : 2.7247e+00 J
  → Total energy    : 9.2279e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=818.62 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.918e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.213e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.173e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.024e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=818.62 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.959e-04
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.058e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.520e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=818.62 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 9.795e-06
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.610e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.978e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=818.62 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.898e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.077e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.010e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.615e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.4170e+02 J
  → Fracture energy : 2.7515e+00 J
  → Total energy    : 9.4445e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=816.38 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.785e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.170e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.085e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.749e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=816.38 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.893e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.059e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.892e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.173e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=816.38 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 9.464e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.569e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.153e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.735e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=816.38 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.732e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.060e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.723e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.447e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.6337e+02 J
  → Fracture energy : 2.7769e+00 J
  → Total energy    : 9.6614e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=814.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.662e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.149e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.023e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.049e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=814.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.831e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.040e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.419e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.565e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=814.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 9.154e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.550e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.845e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.277e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=814.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.577e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.045e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.529e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.187e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.8529e+02 J
  → Fracture energy : 2.8012e+00 J
  → Total energy    : 9.8809e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=812.11 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.546e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.144e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.763e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.666e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=812.11 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.773e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.025e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.078e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.243e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=812.11 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.865e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.533e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.620e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.106e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=812.11 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.433e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.032e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.383e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.065e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0077e+03 J
  → Fracture energy : 2.8246e+00 J
  → Total energy    : 1.0106e+03 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=810.07 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.438e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.149e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.422e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.107e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=810.07 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.719e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.014e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.839e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.743e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=810.07 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.594e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.518e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.461e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.727e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=810.07 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.297e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.020e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.279e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.845e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0308e+03 J
  → Fracture energy : 2.8472e+00 J
  → Total energy    : 1.0337e+03 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=808.10 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.336e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.153e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.193e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.903e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=808.10 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.668e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.003e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.699e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.594e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=808.10 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.340e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.504e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.373e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.653e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=808.10 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.170e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.008e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.221e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.776e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0544e+03 J
  → Fracture energy : 2.8692e+00 J
  → Total energy    : 1.0572e+03 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=806.18 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.240e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.165e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.117e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.503e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=806.18 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.620e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.995e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.640e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.227e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=806.18 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.100e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.491e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.334e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.419e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=806.18 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.050e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.977e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.192e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.651e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.0781e+03 J
  → Fracture energy : 2.8906e+00 J
  → Total energy    : 1.0810e+03 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=804.31 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.150e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.184e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.076e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.188e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=804.31 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.575e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.987e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.626e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.963e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=804.31 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.874e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.479e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.330e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.201e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=804.31 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.937e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.870e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.188e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.482e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1021e+03 J
  → Fracture energy : 2.9119e+00 J
  → Total energy    : 1.1050e+03 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=802.49 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.064e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.198e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.135e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.959e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=802.49 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.532e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.972e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.677e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.767e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=802.49 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.661e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.462e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.366e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.093e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=802.49 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.831e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.745e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.208e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.438e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1265e+03 J
  → Fracture energy : 2.9335e+00 J
  → Total energy    : 1.1295e+03 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=800.72 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.984e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.212e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.193e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.592e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=800.72 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.492e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.956e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.734e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.411e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=800.72 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.459e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.445e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.406e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.848e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=800.72 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.730e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.617e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.233e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.292e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1517e+03 J
  → Fracture energy : 2.9553e+00 J
  → Total energy    : 1.1546e+03 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=799.00 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.907e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.220e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.133e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.430e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=799.00 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.454e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.938e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.702e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.298e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=799.00 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.268e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.427e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.398e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.751e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=799.00 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.634e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.488e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.233e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.210e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.1779e+03 J
  → Fracture energy : 2.9772e+00 J
  → Total energy    : 1.1809e+03 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=797.32 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.835e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.213e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.008e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.355e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=797.32 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.417e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.917e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.611e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.179e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=797.32 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.087e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.409e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.337e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.645e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=797.32 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.543e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.360e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.196e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.141e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2053e+03 J
  → Fracture energy : 2.9990e+00 J
  → Total energy    : 1.2083e+03 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=795.68 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.766e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.186e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.896e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.557e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=795.68 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.383e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.892e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.483e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.423e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=795.68 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.914e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.390e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.240e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.855e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=795.68 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.457e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.231e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.132e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.288e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2332e+03 J
  → Fracture energy : 3.0206e+00 J
  → Total energy    : 1.2363e+03 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=794.08 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.700e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.151e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.652e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.637e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=794.08 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.350e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.864e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.266e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.471e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=794.08 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.750e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.370e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.090e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.893e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=794.08 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.375e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.101e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.037e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.320e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2611e+03 J
  → Fracture energy : 3.0417e+00 J
  → Total energy    : 1.2641e+03 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=792.52 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.637e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.121e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.356e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.502e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=792.52 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.319e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.839e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.001e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.393e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=792.52 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.594e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.351e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.901e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.864e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=792.52 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.297e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.978e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.915e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.316e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.2887e+03 J
  → Fracture energy : 3.0624e+00 J
  → Total energy    : 1.2917e+03 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=790.99 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.578e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.094e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.060e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.568e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=790.99 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.289e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.814e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.740e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.403e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=790.99 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.444e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.333e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.714e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.809e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=790.99 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.222e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.859e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.790e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.284e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3162e+03 J
  → Fracture energy : 3.0827e+00 J
  → Total energy    : 1.3193e+03 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=789.50 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.521e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.061e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.792e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.547e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=789.50 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.260e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.789e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.470e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.371e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=789.50 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.302e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.316e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.509e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.793e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=789.50 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.151e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.742e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.654e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.234e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3436e+03 J
  → Fracture energy : 3.1024e+00 J
  → Total energy    : 1.3467e+03 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=788.04 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.466e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.024e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.467e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.333e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=788.04 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.233e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.763e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.170e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.222e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=788.04 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.166e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.298e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.292e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.713e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=788.04 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.083e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.627e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.513e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.196e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3708e+03 J
  → Fracture energy : 3.1214e+00 J
  → Total energy    : 1.3739e+03 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=786.61 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.414e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.980e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.148e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.217e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=786.61 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.207e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.736e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.877e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.116e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=786.61 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.035e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.280e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.082e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.643e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=786.61 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.018e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.514e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.376e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.158e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.3977e+03 J
  → Fracture energy : 3.1397e+00 J
  → Total energy    : 1.4009e+03 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=785.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.364e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.935e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.872e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.990e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=785.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.182e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.710e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.611e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.877e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=785.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.910e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.263e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.885e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.465e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=785.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.955e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.404e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.243e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.043e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4243e+03 J
  → Fracture energy : 3.1576e+00 J
  → Total energy    : 1.4275e+03 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=783.85 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.316e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.897e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.577e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.772e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=783.85 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.158e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.685e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.334e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.638e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=783.85 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.791e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.246e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.684e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.241e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=783.85 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.895e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.298e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.112e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.906e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4506e+03 J
  → Fracture energy : 3.1750e+00 J
  → Total energy    : 1.4538e+03 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=782.50 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.270e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.858e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.285e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.522e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=782.50 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.135e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.661e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.064e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.485e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=782.50 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.676e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.230e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.486e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.179e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=782.50 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.838e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.197e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.984e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.848e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4766e+03 J
  → Fracture energy : 3.1919e+00 J
  → Total energy    : 1.4798e+03 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=781.19 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.226e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.814e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.016e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.343e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=781.19 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.113e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.636e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.794e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.310e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=781.19 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.565e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.214e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.286e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.056e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=781.19 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.783e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.098e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.851e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.773e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5021e+03 J
  → Fracture energy : 3.2081e+00 J
  → Total energy    : 1.5053e+03 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=779.90 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.184e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.772e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.744e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.087e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=779.90 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.092e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.612e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.539e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.054e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=779.90 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.459e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.199e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.100e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.869e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=779.90 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.730e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.001e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.727e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.651e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5272e+03 J
  → Fracture energy : 3.2239e+00 J
  → Total energy    : 1.5304e+03 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=778.63 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.143e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.734e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.489e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.831e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=778.63 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.071e-04
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.298e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.785e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=778.63 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.357e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.185e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.922e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.687e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=778.63 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.679e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.910e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.609e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.538e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5518e+03 J
  → Fracture energy : 3.2391e+00 J
  → Total energy    : 1.5551e+03 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=777.39 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.104e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.694e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.241e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.644e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=777.39 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.052e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.568e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.056e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.598e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=777.39 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.259e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.171e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.742e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.533e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=777.39 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.629e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.821e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.491e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.440e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5760e+03 J
  → Fracture energy : 3.2538e+00 J
  → Total energy    : 1.5793e+03 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=776.17 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.066e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.654e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.998e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.531e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=776.17 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.033e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.546e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.832e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.498e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=776.17 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.164e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.157e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.579e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.407e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=776.17 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.582e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.736e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.384e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.308e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.5997e+03 J
  → Fracture energy : 3.2681e+00 J
  → Total energy    : 1.6029e+03 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=774.98 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.029e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.618e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.846e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.396e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=774.98 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.015e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.526e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.672e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.371e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=774.98 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.073e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.144e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.455e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.320e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=774.98 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.536e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.654e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.298e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.253e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6229e+03 J
  → Fracture energy : 3.2821e+00 J
  → Total energy    : 1.6262e+03 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=773.80 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.994e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.594e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.703e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.280e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=773.80 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 9.969e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.509e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.529e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.260e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=773.80 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.985e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.132e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.347e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.242e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=773.80 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.492e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.575e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.226e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.206e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6459e+03 J
  → Fracture energy : 3.2958e+00 J
  → Total energy    : 1.6492e+03 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=772.65 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.960e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.579e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.555e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.155e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=772.65 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 9.799e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.493e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.392e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.135e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=772.65 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.899e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.120e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.248e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.152e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=772.65 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.450e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.499e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.162e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.148e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6688e+03 J
  → Fracture energy : 3.3093e+00 J
  → Total energy    : 1.6721e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=771.51 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.927e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.569e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.428e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.983e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=771.51 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 9.635e-05
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.267e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.960e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=771.51 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.817e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.110e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.155e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.021e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=771.51 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.409e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.427e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.100e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.062e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.6915e+03 J
  → Fracture energy : 3.3224e+00 J
  → Total energy    : 1.6948e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=770.40 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.895e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.560e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.299e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.920e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=770.40 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 9.476e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.466e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.153e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.972e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=770.40 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.738e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.099e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.071e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.072e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=770.40 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.369e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.356e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.043e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.121e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7140e+03 J
  → Fracture energy : 3.3353e+00 J
  → Total energy    : 1.7173e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=769.30 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.864e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.550e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.233e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.995e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=769.30 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 9.322e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.453e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.073e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.042e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=769.30 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.661e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.089e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.005e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.127e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=769.30 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.330e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.288e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.995e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.161e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7362e+03 J
  → Fracture energy : 3.3481e+00 J
  → Total energy    : 1.7396e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=768.22 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.835e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.536e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.151e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.999e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=768.22 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 9.173e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.440e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.992e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.037e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=768.22 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.586e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.079e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.943e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.123e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=768.22 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.293e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.219e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.953e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.159e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7582e+03 J
  → Fracture energy : 3.3608e+00 J
  → Total energy    : 1.7615e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=767.16 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.806e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.515e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.040e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.996e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=767.16 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 9.029e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.425e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.885e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.043e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=767.16 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.514e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.069e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.867e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.116e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=767.16 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.257e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.153e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.903e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.143e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7797e+03 J
  → Fracture energy : 3.3732e+00 J
  → Total energy    : 1.7831e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=766.11 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.778e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.488e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.886e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.010e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=766.11 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 8.889e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.409e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.743e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.062e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=766.11 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.445e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.058e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.768e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.139e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=766.11 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.222e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.088e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.840e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.164e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8008e+03 J
  → Fracture energy : 3.3852e+00 J
  → Total energy    : 1.8042e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=765.08 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.751e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.462e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.717e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.965e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=765.08 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 8.754e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.395e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.585e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.020e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=765.08 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.377e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.049e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.656e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.113e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=765.08 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.188e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.027e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.770e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.151e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8214e+03 J
  → Fracture energy : 3.3968e+00 J
  → Total energy    : 1.8248e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=764.07 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.725e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.437e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.582e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.859e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=764.07 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 8.623e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.381e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.473e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.917e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=764.07 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.311e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.040e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.578e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.039e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=764.07 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.156e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.970e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.719e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.105e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8414e+03 J
  → Fracture energy : 3.4082e+00 J
  → Total energy    : 1.8448e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=763.07 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.699e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.415e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.495e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.709e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=763.07 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 8.495e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.368e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.382e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.768e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=763.07 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.248e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.031e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.509e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.928e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=763.07 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.124e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.913e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.673e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.032e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8610e+03 J
  → Fracture energy : 3.4194e+00 J
  → Total energy    : 1.8644e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=762.09 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.674e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.394e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.390e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.636e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=762.09 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 8.372e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.355e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.290e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.701e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=762.09 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.186e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.022e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.446e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.864e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=762.09 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.093e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.857e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.632e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.976e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8801e+03 J
  → Fracture energy : 3.4303e+00 J
  → Total energy    : 1.8835e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=761.12 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.650e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.375e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.310e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.591e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=761.12 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 8.252e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.343e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.209e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.664e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=761.12 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.126e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.014e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.383e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.843e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=761.12 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.063e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.803e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.591e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.967e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.8988e+03 J
  → Fracture energy : 3.4412e+00 J
  → Total energy    : 1.9023e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=760.17 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.627e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.360e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.220e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.524e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=760.17 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 8.135e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.331e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.122e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.602e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=760.17 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.067e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.006e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.320e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.802e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=760.17 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.034e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.750e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.549e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.943e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9173e+03 J
  → Fracture energy : 3.4517e+00 J
  → Total energy    : 1.9208e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=759.23 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.604e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.343e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.112e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.424e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=759.23 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 8.021e-05
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.024e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.506e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=759.23 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.011e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.983e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.250e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.734e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=759.23 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.005e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.700e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.504e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.900e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9355e+03 J
  → Fracture energy : 3.4620e+00 J
  → Total energy    : 1.9390e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=758.30 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.582e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.327e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.020e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.286e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=758.30 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 7.911e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.310e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.947e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.370e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=758.30 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.956e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.911e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.196e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.634e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=758.30 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.978e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.653e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.469e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.836e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9533e+03 J
  → Fracture energy : 3.4722e+00 J
  → Total energy    : 1.9568e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=757.39 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.561e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.311e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.973e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.200e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=757.39 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 7.804e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.300e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.890e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.274e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=757.39 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.902e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.843e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.149e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.538e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=757.39 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.951e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.607e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.437e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.756e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9707e+03 J
  → Fracture energy : 3.4823e+00 J
  → Total energy    : 1.9742e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=756.49 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.540e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.298e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.925e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.181e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=756.49 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 7.700e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.291e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.838e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.263e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=756.49 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.850e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.776e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.110e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.537e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=756.49 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.925e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.563e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.410e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.757e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9877e+03 J
  → Fracture energy : 3.4923e+00 J
  → Total energy    : 1.9912e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=755.60 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.520e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.287e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.885e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.156e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=755.60 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 7.598e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.282e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.801e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.243e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=755.60 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.799e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.712e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.082e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.528e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=755.60 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.899e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.521e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.391e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.756e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0043e+03 J
  → Fracture energy : 3.5023e+00 J
  → Total energy    : 2.0078e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=754.72 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.500e-03
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.857e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.124e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=754.72 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 7.499e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.274e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.768e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.215e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=754.72 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.749e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.649e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.056e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.512e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=754.72 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.875e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.478e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.373e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.748e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0206e+03 J
  → Fracture energy : 3.5122e+00 J
  → Total energy    : 2.0241e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=753.86 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.480e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.271e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.806e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.080e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=753.86 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 7.402e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.266e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.724e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.174e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=753.86 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.701e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.588e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.025e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.484e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=753.86 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.851e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.437e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.353e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.732e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0365e+03 J
  → Fracture energy : 3.5219e+00 J
  → Total energy    : 2.0401e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=753.01 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.462e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.265e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.746e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.019e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=753.01 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 7.309e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.258e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.669e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.114e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=753.01 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.654e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.530e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.986e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.441e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=753.01 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.827e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.397e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.329e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.704e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0522e+03 J
  → Fracture energy : 3.5315e+00 J
  → Total energy    : 2.0558e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=752.16 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.443e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.258e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.694e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.937e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=752.16 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 7.217e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.251e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.628e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.031e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=752.16 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.608e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.474e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.959e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.380e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=752.16 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.804e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.359e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.312e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.664e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0677e+03 J
  → Fracture energy : 3.5409e+00 J
  → Total energy    : 2.0712e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=751.33 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.426e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.251e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.675e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.936e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=751.33 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 7.128e-05
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.607e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.019e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=751.33 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.564e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.418e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.940e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.350e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=751.33 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.782e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.322e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.297e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.628e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0828e+03 J
  → Fracture energy : 3.5503e+00 J
  → Total energy    : 2.0864e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=750.51 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.408e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.244e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.651e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.954e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=750.51 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 7.040e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.237e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.580e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.043e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=750.51 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.520e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.363e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.920e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.374e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=750.51 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.760e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.284e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.283e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.648e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.0978e+03 J
  → Fracture energy : 3.5597e+00 J
  → Total energy    : 2.1013e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=749.70 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.391e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.236e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.610e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.960e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=749.70 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.955e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.229e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.547e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.053e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=749.70 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.478e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.308e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.897e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.386e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=749.70 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.739e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.247e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.268e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.660e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1125e+03 J
  → Fracture energy : 3.5690e+00 J
  → Total energy    : 2.1161e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=748.90 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.374e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.228e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.579e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.955e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=748.90 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.872e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.222e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.512e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.050e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=748.90 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.436e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.254e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.870e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.387e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=748.90 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.718e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.211e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.251e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.664e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1270e+03 J
  → Fracture energy : 3.5782e+00 J
  → Total energy    : 2.1306e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=748.11 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.358e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.220e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.530e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.943e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=748.11 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.791e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.215e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.466e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.040e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=748.11 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.396e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.202e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.838e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.382e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=748.11 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.698e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.176e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.231e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.662e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1414e+03 J
  → Fracture energy : 3.5871e+00 J
  → Total energy    : 2.1450e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=747.33 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.342e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.212e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.473e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.921e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=747.33 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.712e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.209e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.415e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.017e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=747.33 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.356e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.153e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.802e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.366e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=747.33 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.678e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.143e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.208e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.653e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1555e+03 J
  → Fracture energy : 3.5959e+00 J
  → Total energy    : 2.1591e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=746.56 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.327e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.206e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.421e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.877e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=746.56 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.635e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.203e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.372e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.971e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=746.56 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.317e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.107e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.772e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.332e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=746.56 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.659e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.112e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.189e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.630e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1694e+03 J
  → Fracture energy : 3.6046e+00 J
  → Total energy    : 2.1730e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=745.79 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.312e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.200e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.410e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.830e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=745.79 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.559e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.197e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.358e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.910e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=745.79 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.280e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.062e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.760e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.279e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=745.79 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.640e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.081e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.179e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.594e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1831e+03 J
  → Fracture energy : 3.6133e+00 J
  → Total energy    : 2.1867e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=745.04 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.297e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.195e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.403e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.864e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=745.04 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.485e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.191e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.346e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.949e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=745.04 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.243e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.016e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.749e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.299e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=745.04 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.621e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.050e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.171e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.595e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1967e+03 J
  → Fracture energy : 3.6220e+00 J
  → Total energy    : 2.2003e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=744.29 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.283e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.189e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.369e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.894e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=744.29 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.413e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.185e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.316e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.982e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=744.29 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.207e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.972e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.729e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.328e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=744.29 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.603e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.020e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.158e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.618e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2100e+03 J
  → Fracture energy : 3.6306e+00 J
  → Total energy    : 2.2136e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=743.56 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.268e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.183e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.343e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.908e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=743.56 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.342e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.179e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.293e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.999e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=743.56 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.171e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.927e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.713e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.344e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=743.56 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.586e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.990e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.148e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.631e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2232e+03 J
  → Fracture energy : 3.6391e+00 J
  → Total energy    : 2.2268e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=742.83 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.255e-03
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.312e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.910e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=742.83 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.273e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.173e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.265e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.002e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=742.83 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.137e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.882e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.693e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.348e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=742.83 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.568e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.959e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.136e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.636e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2362e+03 J
  → Fracture energy : 3.6476e+00 J
  → Total energy    : 2.2398e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=742.11 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.241e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.168e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.279e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.895e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=742.11 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.206e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.167e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.232e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.988e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=742.11 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.103e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.835e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.670e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.339e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=742.11 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.551e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.929e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.121e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.631e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2490e+03 J
  → Fracture energy : 3.6559e+00 J
  → Total energy    : 2.2527e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=741.40 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.228e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.159e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.235e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.857e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=741.40 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.140e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.160e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.193e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.949e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=741.40 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.070e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.787e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.643e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.310e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=741.40 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.535e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.897e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.104e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.612e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2618e+03 J
  → Fracture energy : 3.6640e+00 J
  → Total energy    : 2.2655e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=740.69 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.215e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.150e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.196e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.791e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=740.69 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.075e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.153e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.166e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.880e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=740.69 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.037e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.739e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.626e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.259e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=740.69 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.519e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.866e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.093e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.578e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2745e+03 J
  → Fracture energy : 3.6721e+00 J
  → Total energy    : 2.2782e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=740.00 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.202e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.140e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.202e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.732e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=740.00 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.011e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.146e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.162e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.814e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=740.00 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.006e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.690e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.619e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.193e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=740.00 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.503e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.835e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.087e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.529e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2871e+03 J
  → Fracture energy : 3.6804e+00 J
  → Total energy    : 2.2908e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=739.31 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.190e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.131e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.191e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.747e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=739.31 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 5.949e-05
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.150e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.833e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=739.31 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.975e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.642e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.610e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.211e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=739.31 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.487e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.804e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.082e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.535e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2996e+03 J
  → Fracture energy : 3.6886e+00 J
  → Total energy    : 2.3033e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=738.63 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.178e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.122e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.175e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.755e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=738.63 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 5.888e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.132e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.137e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.844e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=738.63 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.944e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.596e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.601e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.222e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=738.63 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.472e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.774e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.076e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.546e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3120e+03 J
  → Fracture energy : 3.6968e+00 J
  → Total energy    : 2.3157e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=737.95 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.166e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.114e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.172e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.750e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=737.95 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 5.829e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.126e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.129e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.841e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=737.95 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.914e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.551e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.593e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.223e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=737.95 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.457e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.745e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.070e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.548e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3242e+03 J
  → Fracture energy : 3.7050e+00 J
  → Total energy    : 2.3279e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=737.28 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.154e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.108e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.155e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.734e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=737.28 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 5.770e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.120e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.113e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.827e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=737.28 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.885e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.509e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.582e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.214e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=737.28 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.443e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.717e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.063e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.544e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3364e+03 J
  → Fracture energy : 3.7132e+00 J
  → Total energy    : 2.3401e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=736.62 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.143e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.101e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.127e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.704e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=736.62 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 5.713e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.115e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.089e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.797e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=736.62 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.856e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.470e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.566e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.194e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=736.62 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.428e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.691e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.052e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.531e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3484e+03 J
  → Fracture energy : 3.7212e+00 J
  → Total energy    : 2.3521e+03 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=735.97 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.131e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.094e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.103e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.654e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=735.97 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 5.657e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.110e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.067e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.748e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=735.97 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.828e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.431e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.550e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.157e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=735.97 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.414e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.665e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.042e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.507e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3603e+03 J
  → Fracture energy : 3.7292e+00 J
  → Total energy    : 2.3640e+03 J
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
  T_new: min=263.15 K, max=1023.18 K, mean=735.32 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.120e-03
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.096e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.584e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=735.32 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 5.602e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.104e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.059e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.676e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=735.32 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.801e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.391e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.542e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.102e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=735.32 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.400e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.639e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.036e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.470e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3720e+03 J
  → Fracture energy : 3.7372e+00 J
  → Total energy    : 2.3757e+03 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0099.vtu

Simulation completed in 311.31 s
Total time steps solved: 100
