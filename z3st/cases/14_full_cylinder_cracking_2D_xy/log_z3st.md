Info    : Reading 'mesh.msh'...
Info    : 17968 nodes
Info    : 35954 elements
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
  split               : star_convex
  gamma_star          : 0.0
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
    Set 17968 DOFs to 1023.15 K
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
Computing symbolic result fields (strain, stress, ...)

**[INFO]** Hot-reload of allow-listed input.yaml parameters is active. Edit input.yaml during the run; changes apply at the next step boundary. Allowed keys: damage.{stag_tol,rtol,hybrid_constraint,gamma_star}, mechanical.{stag_tol,rtol}, thermal.{stag_tol,rtol}, solver_settings.{max_iters,relax_*}.


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
  T_new: min=263.15 K, max=1023.15 K, mean=996.33 K
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
  T_new: min=263.15 K, max=1023.15 K, mean=996.33 K
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
  T_new: min=263.15 K, max=1023.15 K, mean=996.33 K
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
  T_new: min=263.15 K, max=1023.15 K, mean=996.33 K
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
  T_new: min=263.15 K, max=1023.15 K, mean=996.33 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.127e-05
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
  T_new: min=263.15 K, max=1023.15 K, mean=996.33 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.637e-07
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
  T_new: min=263.15 K, max=1023.15 K, mean=996.33 K
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
  T_new: min=263.15 K, max=1023.15 K, mean=996.33 K
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
  T_new: min=263.15 K, max=1023.15 K, mean=996.33 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.046e-11
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
  T_new: min=263.15 K, max=1023.15 K, mean=996.33 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.523e-12
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
  T_new: min=263.15 K, max=1023.15 K, mean=996.33 K
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
  T_new: min=263.15 K, max=1023.15 K, mean=996.33 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 8.808e-15
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
  T_new: min=263.15 K, max=1023.15 K, mean=996.33 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.405e-16
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.024e-12
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
  → Elastic energy  : 2.3519e+02 J
  → Fracture energy : 3.1642e-01 J
  → Total energy    : 2.3551e+02 J


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
  T_new: min=263.15 K, max=1023.15 K, mean=959.69 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.207e-02
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
  T_new: min=263.15 K, max=1023.15 K, mean=959.69 K
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
  T_new: min=263.15 K, max=1023.15 K, mean=959.69 K
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
  T_new: min=263.15 K, max=1023.15 K, mean=959.69 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.151e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.278e-07
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
  T_new: min=263.15 K, max=1023.15 K, mean=959.69 K
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
  → Elastic energy  : 2.3455e+02 J
  → Fracture energy : 7.6716e-01 J
  → Total energy    : 2.3531e+02 J


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
  T_new: min=263.15 K, max=1023.15 K, mean=939.44 K
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
  T_new: min=263.15 K, max=1023.15 K, mean=939.44 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.145e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.446e-05
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
  T_new: min=263.15 K, max=1023.15 K, mean=939.44 K
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
  T_new: min=263.15 K, max=1023.15 K, mean=939.44 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.361e-06
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
  T_new: min=263.15 K, max=1023.15 K, mean=939.44 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.681e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.368e-08
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
  → Elastic energy  : 2.3453e+02 J
  → Fracture energy : 1.0278e+00 J
  → Total energy    : 2.3556e+02 J


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
  T_new: min=263.15 K, max=1023.15 K, mean=925.07 K
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
  T_new: min=263.15 K, max=1023.15 K, mean=925.07 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.424e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.361e-05
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
  T_new: min=263.15 K, max=1023.15 K, mean=925.07 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.119e-05
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
  T_new: min=263.15 K, max=1023.15 K, mean=925.07 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.560e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.017e-07
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
  T_new: min=263.15 K, max=1023.15 K, mean=925.07 K
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
  → Elastic energy  : 2.3443e+02 J
  → Fracture energy : 1.2138e+00 J
  → Total energy    : 2.3564e+02 J


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
  T_new: min=263.15 K, max=1023.15 K, mean=913.74 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.158e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.273e-04
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
  T_new: min=263.15 K, max=1023.15 K, mean=913.74 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.079e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.880e-05
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
  T_new: min=263.15 K, max=1023.15 K, mean=913.74 K
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
  T_new: min=263.15 K, max=1023.15 K, mean=913.74 K
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
  T_new: min=263.15 K, max=1023.15 K, mean=913.74 K
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
  → Elastic energy  : 2.3432e+02 J
  → Fracture energy : 1.3646e+00 J
  → Total energy    : 2.3568e+02 J


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
  T_new: min=263.15 K, max=1023.15 K, mean=904.32 K
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
  T_new: min=263.15 K, max=1023.15 K, mean=904.32 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 8.753e-04
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
  T_new: min=263.15 K, max=1023.15 K, mean=904.32 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.376e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.569e-06
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
  T_new: min=263.15 K, max=1023.15 K, mean=904.32 K
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
  T_new: min=263.15 K, max=1023.15 K, mean=904.32 K
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
  → Elastic energy  : 2.3421e+02 J
  → Fracture energy : 1.4944e+00 J
  → Total energy    : 2.3571e+02 J


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
  T_new: min=263.15 K, max=1023.15 K, mean=896.21 K
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
  T_new: min=263.15 K, max=1023.15 K, mean=896.21 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.396e-04
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
  T_new: min=263.15 K, max=1023.15 K, mean=896.21 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.698e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.349e-06
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
  T_new: min=263.15 K, max=1023.15 K, mean=896.21 K
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
  T_new: min=263.15 K, max=1023.15 K, mean=896.21 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.245e-08
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
  → Elastic energy  : 2.3411e+02 J
  → Fracture energy : 1.6090e+00 J
  → Total energy    : 2.3572e+02 J


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
  T_new: min=263.15 K, max=1023.15 K, mean=889.07 K
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
  T_new: min=263.15 K, max=1023.15 K, mean=889.07 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.423e-04
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
  T_new: min=263.15 K, max=1023.15 K, mean=889.07 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.211e-05
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
  T_new: min=263.15 K, max=1023.15 K, mean=889.07 K
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
  T_new: min=263.15 K, max=1023.15 K, mean=889.07 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 8.029e-08
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
  → Elastic energy  : 2.3401e+02 J
  → Fracture energy : 1.7113e+00 J
  → Total energy    : 2.3573e+02 J


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
  T_new: min=263.15 K, max=1023.15 K, mean=882.67 K
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
  T_new: min=263.15 K, max=1023.15 K, mean=882.67 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.688e-04
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
  T_new: min=263.15 K, max=1023.15 K, mean=882.67 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.844e-05
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
  T_new: min=263.15 K, max=1023.15 K, mean=882.67 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.422e-06
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
  T_new: min=263.15 K, max=1023.15 K, mean=882.67 K
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
  → Elastic energy  : 2.3392e+02 J
  → Fracture energy : 1.8041e+00 J
  → Total energy    : 2.3573e+02 J


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
  T_new: min=263.15 K, max=1023.15 K, mean=876.87 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.023e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.771e-04
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
  T_new: min=263.15 K, max=1023.15 K, mean=876.87 K
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
  T_new: min=263.15 K, max=1023.15 K, mean=876.87 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.556e-05
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
  T_new: min=263.15 K, max=1023.15 K, mean=876.87 K
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
  T_new: min=263.15 K, max=1023.15 K, mean=876.87 K
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
  → Elastic energy  : 2.3383e+02 J
  → Fracture energy : 1.8904e+00 J
  → Total energy    : 2.3572e+02 J


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
  T_new: min=263.15 K, max=1023.15 K, mean=871.56 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.296e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.012e-04
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
  T_new: min=263.15 K, max=1023.15 K, mean=871.56 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.648e-04
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
  T_new: min=263.15 K, max=1023.15 K, mean=871.56 K
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
  |ΔD|_∞ = 9.989e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=871.56 K
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
  T_new: min=263.15 K, max=1023.15 K, mean=871.56 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.810e-08
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
  → Elastic energy  : 2.3374e+02 J
  → Fracture energy : 1.9711e+00 J
  → Total energy    : 2.3571e+02 J


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
  T_new: min=263.15 K, max=1023.15 K, mean=866.66 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 8.530e-03
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
  T_new: min=263.15 K, max=1023.15 K, mean=866.66 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.265e-04
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
  T_new: min=263.15 K, max=1023.15 K, mean=866.66 K
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
  T_new: min=263.15 K, max=1023.15 K, mean=866.66 K
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
  T_new: min=263.15 K, max=1023.15 K, mean=866.66 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.331e-08
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
  → Elastic energy  : 2.3364e+02 J
  → Fracture energy : 2.0484e+00 J
  → Total energy    : 2.3569e+02 J


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
  T_new: min=263.15 K, max=1023.15 K, mean=862.11 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.886e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.969e-04
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
  T_new: min=263.15 K, max=1023.15 K, mean=862.11 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.943e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.825e-05
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
  T_new: min=263.15 K, max=1023.15 K, mean=862.11 K
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
  T_new: min=263.15 K, max=1023.15 K, mean=862.11 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.858e-07
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
  → Elastic energy  : 2.3355e+02 J
  → Fracture energy : 2.1239e+00 J
  → Total energy    : 2.3567e+02 J


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
  T_new: min=263.15 K, max=1023.15 K, mean=857.86 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.337e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.844e-04
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
  T_new: min=263.15 K, max=1023.15 K, mean=857.86 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.669e-04
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
  T_new: min=263.15 K, max=1023.15 K, mean=857.86 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.834e-05
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
  T_new: min=263.15 K, max=1023.15 K, mean=857.86 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.172e-07
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
  → Elastic energy  : 2.3345e+02 J
  → Fracture energy : 2.1953e+00 J
  → Total energy    : 2.3564e+02 J


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
  T_new: min=263.15 K, max=1023.15 K, mean=853.87 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.863e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.808e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.461e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.722e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=853.87 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.431e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.483e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.994e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.525e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=853.87 K
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
  T_new: min=263.15 K, max=1023.15 K, mean=853.87 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 8.579e-07
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
  ||ΔD||/||D|| = 1.385e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.925e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3336e+02 J
  → Fracture energy : 2.2615e+00 J
  → Total energy    : 2.3562e+02 J


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
  T_new: min=263.15 K, max=1023.15 K, mean=850.10 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.449e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.377e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.233e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.796e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=850.10 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.224e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.228e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.809e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.593e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=850.10 K
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
  ||ΔD||/||D|| = 1.978e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.138e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=850.10 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 8.061e-07
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
  |ΔD|_∞ = 7.400e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3327e+02 J
  → Fracture energy : 2.3222e+00 J
  → Total energy    : 2.3559e+02 J


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
  T_new: min=263.15 K, max=1023.15 K, mean=846.54 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.084e-03
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
  |ΔD|_∞ = 1.815e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=846.54 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.042e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.020e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.635e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.635e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=846.54 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.521e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.194e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.853e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.175e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=846.54 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.605e-07
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
  |ΔD|_∞ = 7.658e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3318e+02 J
  → Fracture energy : 2.3787e+00 J
  → Total energy    : 2.3556e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=843.17 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.760e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.732e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.799e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.767e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=843.17 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.880e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.845e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.450e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.602e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=843.17 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.440e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.085e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.725e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.163e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=843.17 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.200e-07
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
  |ΔD|_∞ = 7.662e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3310e+02 J
  → Fracture energy : 2.4309e+00 J
  → Total energy    : 2.3553e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=839.95 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.470e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.517e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.572e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.648e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=839.95 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.735e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.709e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.258e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.508e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=839.95 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.367e-05
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
  |ΔD|_∞ = 1.105e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=839.95 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.837e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.343e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.031e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.347e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3302e+02 J
  → Fracture energy : 2.4789e+00 J
  → Total energy    : 2.3550e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=836.89 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.209e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.350e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.339e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.567e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=836.89 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.604e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.604e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.061e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.394e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=836.89 K
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
  |ΔD|_∞ = 1.026e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=836.89 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.511e-07
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
  ||ΔD||/||D|| = 9.405e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.879e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3295e+02 J
  → Fracture energy : 2.5229e+00 J
  → Total energy    : 2.3547e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=833.96 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.972e-03
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
  ||ΔD||/||D|| = 2.122e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.513e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=833.96 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.486e-04
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
  |ΔD|_∞ = 1.356e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=833.96 K
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
  |ΔD|_∞ = 9.606e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=833.96 K
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
  ||ΔD||/||D|| = 8.548e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.271e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3288e+02 J
  → Fracture energy : 2.5632e+00 J
  → Total energy    : 2.3544e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=831.15 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.757e-03
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
  |ΔD|_∞ = 1.438e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=831.15 K
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
  ||ΔD||/||D|| = 1.692e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.300e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=831.15 K
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
  ||ΔD||/||D|| = 1.193e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.242e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=831.15 K
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
  ||ΔD||/||D|| = 7.746e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.922e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3281e+02 J
  → Fracture energy : 2.6005e+00 J
  → Total energy    : 2.3541e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=828.46 K
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
  T_new: min=263.15 K, max=1023.16 K, mean=828.46 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.280e-04
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
  ||ΔD||/||D|| = 1.531e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.222e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=828.46 K
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
  ||ΔD||/||D|| = 1.081e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.721e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=828.46 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.701e-07
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
  ||ΔD||/||D|| = 7.022e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.662e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3275e+02 J
  → Fracture energy : 2.6350e+00 J
  → Total energy    : 2.3538e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=825.87 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.380e-03
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
  ||ΔD||/||D|| = 1.563e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.269e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=825.87 K
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
  T_new: min=263.15 K, max=1023.16 K, mean=825.87 K
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
  ||ΔD||/||D|| = 9.810e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.309e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=825.87 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.475e-07
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
  ||ΔD||/||D|| = 6.389e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.383e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3268e+02 J
  → Fracture energy : 2.6668e+00 J
  → Total energy    : 2.3535e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=823.38 K
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
  |ΔD|_∞ = 1.200e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=823.38 K
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
  |ΔD|_∞ = 1.105e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=823.38 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.053e-05
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
  ||ΔD||/||D|| = 8.938e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.018e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=823.38 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.267e-07
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
  ||ΔD||/||D|| = 5.841e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.250e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3263e+02 J
  → Fracture energy : 2.6966e+00 J
  → Total energy    : 2.3532e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=820.98 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.060e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.285e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.281e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.111e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=820.98 K
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
  ||ΔD||/||D|| = 1.148e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.032e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=820.98 K
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
  ||ΔD||/||D|| = 8.201e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.482e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=820.98 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.075e-07
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
  ||ΔD||/||D|| = 5.377e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.881e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3257e+02 J
  → Fracture energy : 2.7249e+00 J
  → Total energy    : 2.3529e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=818.66 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.918e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.216e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.175e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.027e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=818.66 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.959e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.086e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.060e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.547e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=818.66 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 9.794e-06
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
  ||ΔD||/||D|| = 7.618e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.996e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=818.66 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.897e-07
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
  ||ΔD||/||D|| = 5.014e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.626e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3251e+02 J
  → Fracture energy : 2.7516e+00 J
  → Total energy    : 2.3526e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=816.42 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.785e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.173e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.086e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.776e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=816.42 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.892e-04
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
  ||ΔD||/||D|| = 9.905e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.196e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=816.42 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 9.462e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.568e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.160e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.749e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=816.42 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.731e-07
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
  ||ΔD||/||D|| = 4.726e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.455e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3246e+02 J
  → Fracture energy : 2.7770e+00 J
  → Total energy    : 2.3524e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=814.25 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.661e-03
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
  ||ΔD||/||D|| = 1.024e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.070e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=814.25 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.831e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.039e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.429e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.582e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=814.25 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 9.153e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.549e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.851e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.293e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=814.25 K
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
  ||ΔD||/||D|| = 4.532e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.198e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3241e+02 J
  → Fracture energy : 2.8013e+00 J
  → Total energy    : 2.3521e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=812.15 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.546e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.137e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.767e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.679e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=812.15 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.773e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.023e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.081e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.255e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=812.15 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.864e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.532e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.621e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.113e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=812.15 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.432e-07
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
  |ΔD|_∞ = 4.069e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3235e+02 J
  → Fracture energy : 2.8248e+00 J
  → Total energy    : 2.3518e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=810.11 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.437e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.142e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.418e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.100e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=810.11 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.719e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.011e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.835e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.735e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=810.11 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.593e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.517e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.458e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.730e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=810.11 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.297e-07
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.277e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.847e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3230e+02 J
  → Fracture energy : 2.8474e+00 J
  → Total energy    : 2.3515e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=808.13 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.335e-03
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
  ||ΔD||/||D|| = 9.185e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.893e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=808.13 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.668e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.002e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.691e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.584e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=808.13 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.339e-06
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
  ||ΔD||/||D|| = 6.368e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.644e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=808.13 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.169e-07
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
  ||ΔD||/||D|| = 4.218e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.770e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3225e+02 J
  → Fracture energy : 2.8693e+00 J
  → Total energy    : 2.3512e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=806.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.240e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.167e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.107e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.494e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=806.21 K
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
  ||ΔD||/||D|| = 8.632e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.220e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=806.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.099e-06
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
  ||ΔD||/||D|| = 6.328e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.413e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=806.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.050e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.976e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.189e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.647e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3220e+02 J
  → Fracture energy : 2.8908e+00 J
  → Total energy    : 2.3509e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=804.35 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.149e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.183e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.066e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.168e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=804.35 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.575e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.986e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.617e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.943e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=804.35 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.873e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.478e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.323e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.186e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=804.35 K
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
  ||ΔD||/||D|| = 4.183e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.471e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3215e+02 J
  → Fracture energy : 2.9120e+00 J
  → Total energy    : 2.3506e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=802.53 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.064e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.192e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.123e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.937e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=802.53 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.532e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.970e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.666e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.747e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=802.53 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.660e-06
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
  ||ΔD||/||D|| = 6.358e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.078e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=802.53 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.830e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.743e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.203e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.427e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3210e+02 J
  → Fracture energy : 2.9335e+00 J
  → Total energy    : 2.3504e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=800.76 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.983e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.199e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.181e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.574e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=800.76 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.492e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.952e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.724e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.397e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=800.76 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.458e-06
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
  ||ΔD||/||D|| = 6.400e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.837e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=800.76 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.729e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.613e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.229e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.284e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3205e+02 J
  → Fracture energy : 2.9554e+00 J
  → Total energy    : 2.3501e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=799.04 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.907e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.199e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.127e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.407e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=799.04 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.453e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.932e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.696e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.275e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=799.04 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.267e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.425e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.394e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.733e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=799.04 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.634e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.480e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.231e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.197e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3200e+02 J
  → Fracture energy : 2.9772e+00 J
  → Total energy    : 2.3498e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=797.36 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.834e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.192e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.009e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.361e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=797.36 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.417e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.910e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.611e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.184e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=797.36 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.086e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.406e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.337e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.649e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=797.36 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.543e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.350e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.196e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.137e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3195e+02 J
  → Fracture energy : 2.9991e+00 J
  → Total energy    : 2.3495e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=795.72 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.765e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.177e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.898e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.563e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=795.72 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.383e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.889e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.486e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.430e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=795.72 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.913e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.388e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.242e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.859e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=795.72 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.457e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.225e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.133e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.291e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3190e+02 J
  → Fracture energy : 3.0207e+00 J
  → Total energy    : 2.3492e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=794.12 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.700e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.157e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.660e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.641e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=794.12 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.350e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.866e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.272e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.474e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=794.12 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.749e-06
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
  ||ΔD||/||D|| = 6.094e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.895e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=794.12 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.375e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.102e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.039e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.321e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3185e+02 J
  → Fracture energy : 3.0419e+00 J
  → Total energy    : 2.3489e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=792.56 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.637e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.139e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.361e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.503e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=792.56 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.319e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.843e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.004e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.396e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=792.56 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.593e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.353e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.904e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.867e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=792.56 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.296e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.983e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.917e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.318e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3180e+02 J
  → Fracture energy : 3.0626e+00 J
  → Total energy    : 2.3487e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=791.03 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.577e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.111e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.065e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.572e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=791.03 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.289e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.819e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.744e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.407e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=791.03 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.444e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.335e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.716e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.811e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=791.03 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.222e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.865e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.791e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.285e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3176e+02 J
  → Fracture energy : 3.0828e+00 J
  → Total energy    : 2.3484e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=789.54 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.520e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.071e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.797e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.549e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=789.54 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.260e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.792e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.473e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.374e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=789.54 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.301e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.317e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.511e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.795e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=789.54 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.151e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.745e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.655e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.235e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3171e+02 J
  → Fracture energy : 3.1025e+00 J
  → Total energy    : 2.3481e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=788.08 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.466e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.029e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.471e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.333e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=788.08 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.233e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.764e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.174e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.225e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=788.08 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.165e-06
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
  ||ΔD||/||D|| = 5.295e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.716e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=788.08 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.082e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.628e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.514e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.198e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3166e+02 J
  → Fracture energy : 3.1215e+00 J
  → Total energy    : 2.3479e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=786.65 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.414e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.982e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.153e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.218e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=786.65 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.207e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.737e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.881e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.117e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=786.65 K
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
  ||ΔD||/||D|| = 5.084e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.644e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=786.65 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.017e-07
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
  ||ΔD||/||D|| = 3.377e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.159e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3162e+02 J
  → Fracture energy : 3.1399e+00 J
  → Total energy    : 2.3476e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=785.25 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.364e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.937e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.873e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.989e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=785.25 K
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
  ||ΔD||/||D|| = 6.613e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.875e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=785.25 K
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
  ||ΔD||/||D|| = 4.886e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.464e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=785.25 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.955e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.403e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.243e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.042e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3157e+02 J
  → Fracture energy : 3.1578e+00 J
  → Total energy    : 2.3473e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=783.89 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.316e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.898e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.580e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.774e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=783.89 K
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
  ||ΔD||/||D|| = 6.337e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.640e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=783.89 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.790e-06
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
  ||ΔD||/||D|| = 4.685e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.241e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=783.89 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.895e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.297e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.113e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.906e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3153e+02 J
  → Fracture energy : 3.1752e+00 J
  → Total energy    : 2.3471e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=782.54 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.270e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.859e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.286e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.524e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=782.54 K
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
  |ΔD|_∞ = 5.486e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=782.54 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.675e-06
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
  ||ΔD||/||D|| = 4.487e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.181e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=782.54 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.838e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.196e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.984e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.850e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3149e+02 J
  → Fracture energy : 3.1920e+00 J
  → Total energy    : 2.3468e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=781.23 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.226e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.815e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.019e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.343e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=781.23 K
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
  ||ΔD||/||D|| = 5.795e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.310e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=781.23 K
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
  |ΔD|_∞ = 4.057e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=781.23 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.782e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.097e-08
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
  → Elastic energy  : 2.3144e+02 J
  → Fracture energy : 3.2083e+00 J
  → Total energy    : 2.3465e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=779.94 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.183e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.773e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.747e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.085e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=779.94 K
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
  ||ΔD||/||D|| = 5.541e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.053e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=779.94 K
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
  ||ΔD||/||D|| = 4.101e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.868e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=779.94 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.729e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.000e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.728e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.651e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3140e+02 J
  → Fracture energy : 3.2241e+00 J
  → Total energy    : 2.3463e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=778.67 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.143e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.735e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.491e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.834e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=778.67 K
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
  ||ΔD||/||D|| = 5.300e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.785e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=778.67 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.356e-06
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
  ||ΔD||/||D|| = 3.923e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.687e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=778.67 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.678e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.909e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.610e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.538e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3136e+02 J
  → Fracture energy : 3.2393e+00 J
  → Total energy    : 2.3460e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=777.43 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.103e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.695e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.241e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.648e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=777.43 K
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
  ||ΔD||/||D|| = 5.057e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.603e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=777.43 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.258e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.170e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.743e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.532e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=777.43 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.629e-07
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.491e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.439e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3132e+02 J
  → Fracture energy : 3.2540e+00 J
  → Total energy    : 2.3458e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=776.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.065e-03
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
  ||ΔD||/||D|| = 4.996e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.535e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=776.21 K
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
  ||ΔD||/||D|| = 4.830e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.502e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=776.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.163e-06
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
  ||ΔD||/||D|| = 3.577e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.410e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=776.21 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.582e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.735e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.382e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.308e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3128e+02 J
  → Fracture energy : 3.2683e+00 J
  → Total energy    : 2.3455e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=775.02 K
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
  ||ΔD||/||D|| = 4.847e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.399e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=775.02 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.014e-04
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
  |ΔD|_∞ = 4.375e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=775.02 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.072e-06
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
  ||ΔD||/||D|| = 3.454e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.324e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=775.02 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.536e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.653e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.297e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.256e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3124e+02 J
  → Fracture energy : 3.2823e+00 J
  → Total energy    : 2.3452e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=773.84 K
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
  ||ΔD||/||D|| = 4.702e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.284e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=773.84 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 9.968e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.508e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.527e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.264e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=773.84 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.984e-06
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
  ||ΔD||/||D|| = 3.345e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.246e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=773.84 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.492e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.574e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.225e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.208e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3120e+02 J
  → Fracture energy : 3.2960e+00 J
  → Total energy    : 2.3450e+02 J


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
  T_new: min=263.15 K, max=1023.16 K, mean=772.69 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.960e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.578e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.550e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.160e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=772.69 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 9.798e-05
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
  ||ΔD||/||D|| = 4.387e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.140e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=772.69 K
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
  ||ΔD||/||D|| = 3.244e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.156e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=772.69 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.449e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.498e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.159e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.151e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3116e+02 J
  → Fracture energy : 3.3094e+00 J
  → Total energy    : 2.3447e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=771.55 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.927e-03
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.422e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.987e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=771.55 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 9.633e-05
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
  ||ΔD||/||D|| = 4.261e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.964e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=771.55 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.817e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.109e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.151e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.024e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=771.55 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.408e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.426e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.097e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.064e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3112e+02 J
  → Fracture energy : 3.3225e+00 J
  → Total energy    : 2.3445e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=770.44 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.895e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.559e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.294e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.897e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=770.44 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 9.474e-05
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
  ||ΔD||/||D|| = 4.148e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.950e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=770.44 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.737e-06
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
  ||ΔD||/||D|| = 3.067e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.055e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=770.44 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.369e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.355e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.041e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.109e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3109e+02 J
  → Fracture energy : 3.3354e+00 J
  → Total energy    : 2.3442e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=769.34 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.864e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.549e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.229e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.980e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=769.34 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 9.321e-05
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
  ||ΔD||/||D|| = 4.069e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.027e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=769.34 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.660e-06
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
  ||ΔD||/||D|| = 3.003e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.116e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=769.34 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.330e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.287e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.994e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.152e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3105e+02 J
  → Fracture energy : 3.3482e+00 J
  → Total energy    : 2.3440e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=768.26 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.834e-03
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.155e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.996e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=768.26 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 9.172e-05
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
  ||ΔD||/||D|| = 3.995e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.034e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=768.26 K
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
  ||ΔD||/||D|| = 2.945e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.121e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=768.26 K
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
  |ΔD|_∞ = 2.158e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3101e+02 J
  → Fracture energy : 3.3609e+00 J
  → Total energy    : 2.3437e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=767.20 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.806e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.517e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.045e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.995e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=767.20 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 9.028e-05
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
  ||ΔD||/||D|| = 3.891e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.041e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=767.20 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.514e-06
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.870e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.114e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=767.20 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.257e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.152e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.905e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.141e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3097e+02 J
  → Fracture energy : 3.3734e+00 J
  → Total energy    : 2.3435e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=766.16 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.778e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.490e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.894e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.017e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=766.16 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 8.888e-05
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
  ||ΔD||/||D|| = 3.750e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.068e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=766.16 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.444e-06
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
  ||ΔD||/||D|| = 2.772e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.142e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=766.16 K
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
  ||ΔD||/||D|| = 1.843e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.165e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3094e+02 J
  → Fracture energy : 3.3854e+00 J
  → Total energy    : 2.3432e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=765.13 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.751e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.463e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.724e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.978e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=765.13 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 8.753e-05
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
  ||ΔD||/||D|| = 3.592e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.032e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=765.13 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.376e-06
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
  ||ΔD||/||D|| = 2.662e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.121e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=765.13 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.188e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.026e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.774e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.156e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3090e+02 J
  → Fracture energy : 3.3970e+00 J
  → Total energy    : 2.3430e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=764.12 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.724e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.438e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.586e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.876e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=764.12 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 8.621e-05
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
  ||ΔD||/||D|| = 3.478e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.932e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=764.12 K
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
  ||ΔD||/||D|| = 2.582e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.050e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=764.12 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.155e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.969e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.721e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.111e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3086e+02 J
  → Fracture energy : 3.4084e+00 J
  → Total energy    : 2.3427e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=763.12 K
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
  ||ΔD||/||D|| = 3.502e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.727e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=763.12 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 8.494e-05
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
  ||ΔD||/||D|| = 3.388e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.784e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=763.12 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.247e-06
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
  ||ΔD||/||D|| = 2.513e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.940e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=763.12 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.124e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.912e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.675e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.040e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3083e+02 J
  → Fracture energy : 3.4196e+00 J
  → Total energy    : 2.3425e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=762.14 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.674e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.395e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.396e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.644e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=762.14 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 8.370e-05
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
  ||ΔD||/||D|| = 3.295e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.708e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=762.14 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.185e-06
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
  ||ΔD||/||D|| = 2.449e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.868e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=762.14 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.093e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.856e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.635e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.978e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3079e+02 J
  → Fracture energy : 3.4306e+00 J
  → Total energy    : 2.3422e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=761.17 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.650e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.376e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.313e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.600e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=761.17 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 8.250e-05
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
  ||ΔD||/||D|| = 3.212e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.672e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=761.17 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.125e-06
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
  ||ΔD||/||D|| = 2.386e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.849e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=761.17 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.063e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.802e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.592e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.970e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3076e+02 J
  → Fracture energy : 3.4414e+00 J
  → Total energy    : 2.3420e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=760.22 K
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
  ||ΔD||/||D|| = 3.223e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.534e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=760.22 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 8.134e-05
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
  ||ΔD||/||D|| = 3.125e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.611e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=760.22 K
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
  ||ΔD||/||D|| = 2.322e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.808e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=760.22 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.033e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.749e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.551e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.947e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3072e+02 J
  → Fracture energy : 3.4520e+00 J
  → Total energy    : 2.3418e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=759.28 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.604e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.344e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.115e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.436e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=759.28 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 8.020e-05
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
  ||ΔD||/||D|| = 3.026e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.516e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=759.28 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.010e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.982e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.252e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.741e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=759.28 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.005e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.699e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.505e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.905e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3069e+02 J
  → Fracture energy : 3.4623e+00 J
  → Total energy    : 2.3415e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=758.35 K
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
  |ΔD|_∞ = 3.298e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=758.35 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 7.910e-05
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
  ||ΔD||/||D|| = 2.945e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.382e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=758.35 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.955e-06
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
  ||ΔD||/||D|| = 2.195e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.643e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=758.35 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.978e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.652e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.469e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.841e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3066e+02 J
  → Fracture energy : 3.4724e+00 J
  → Total energy    : 2.3413e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=757.44 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.561e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.312e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.975e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.205e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=757.44 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 7.803e-05
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
  ||ΔD||/||D|| = 2.891e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.278e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=757.44 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.901e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.842e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.149e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.541e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=757.44 K
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
  |ΔD|_∞ = 1.762e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3062e+02 J
  → Fracture energy : 3.4825e+00 J
  → Total energy    : 2.3411e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=756.54 K
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
  ||ΔD||/||D|| = 2.926e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.187e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=756.54 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 7.698e-05
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
  ||ΔD||/||D|| = 2.839e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.268e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=756.54 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.849e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.775e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.110e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.540e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=756.54 K
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
  ||ΔD||/||D|| = 1.411e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.759e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3059e+02 J
  → Fracture energy : 3.4925e+00 J
  → Total energy    : 2.3408e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=755.65 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.519e-03
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
  ||ΔD||/||D|| = 2.888e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.162e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=755.65 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 7.597e-05
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
  ||ΔD||/||D|| = 2.802e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.248e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=755.65 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.798e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.711e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.082e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.532e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=755.65 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.899e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.520e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.391e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.758e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3056e+02 J
  → Fracture energy : 3.5025e+00 J
  → Total energy    : 2.3406e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=754.77 K
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
  ||ΔD||/||D|| = 2.859e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.130e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=754.77 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 7.498e-05
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
  ||ΔD||/||D|| = 2.770e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.221e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=754.77 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.749e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.648e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.057e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.515e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=754.77 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.874e-07
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
  |ΔD|_∞ = 1.750e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3052e+02 J
  → Fracture energy : 3.5125e+00 J
  → Total energy    : 2.3404e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=753.91 K
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
  ||ΔD||/||D|| = 2.809e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.087e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=753.91 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 7.401e-05
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
  ||ΔD||/||D|| = 2.726e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.180e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=753.91 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.701e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.587e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.026e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.488e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=753.91 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.850e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.436e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.354e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.734e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3049e+02 J
  → Fracture energy : 3.5222e+00 J
  → Total energy    : 2.3401e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=753.05 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.461e-03
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
  ||ΔD||/||D|| = 2.749e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.027e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=753.05 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 7.307e-05
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
  ||ΔD||/||D|| = 2.671e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.121e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=753.05 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.654e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.529e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.988e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.446e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=753.05 K
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
  ||ΔD||/||D|| = 1.330e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.708e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3046e+02 J
  → Fracture energy : 3.5317e+00 J
  → Total energy    : 2.3399e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=752.21 K
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
  ||ΔD||/||D|| = 2.696e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.945e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=752.21 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 7.216e-05
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
  ||ΔD||/||D|| = 2.629e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.039e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=752.21 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.608e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.473e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.959e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.385e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=752.21 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.804e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.358e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.311e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.668e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3043e+02 J
  → Fracture energy : 3.5412e+00 J
  → Total energy    : 2.3397e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=751.38 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.425e-03
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
  ||ΔD||/||D|| = 2.676e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.939e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=751.38 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 7.127e-05
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
  ||ΔD||/||D|| = 2.608e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.021e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=751.38 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.563e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.417e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.941e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.351e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=751.38 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.782e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.321e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.298e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.629e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3039e+02 J
  → Fracture energy : 3.5506e+00 J
  → Total energy    : 2.3394e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=750.56 K
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
  ||ΔD||/||D|| = 2.653e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.958e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=750.56 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 7.039e-05
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
  ||ΔD||/||D|| = 2.582e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.046e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=750.56 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.520e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.362e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.921e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.375e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=750.56 K
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
  ||ΔD||/||D|| = 1.284e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.649e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3036e+02 J
  → Fracture energy : 3.5600e+00 J
  → Total energy    : 2.3392e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=749.75 K
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
  ||ΔD||/||D|| = 2.611e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.964e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=749.75 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.954e-05
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
  ||ΔD||/||D|| = 2.548e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.056e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=749.75 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.477e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.307e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.897e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.388e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=749.75 K
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
  ||ΔD||/||D|| = 1.269e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.661e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3033e+02 J
  → Fracture energy : 3.5693e+00 J
  → Total energy    : 2.3390e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=748.95 K
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
  ||ΔD||/||D|| = 2.581e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.959e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=748.95 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.871e-05
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
  ||ΔD||/||D|| = 2.514e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.054e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=748.95 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.436e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.253e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.871e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.390e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=748.95 K
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
  |ΔD|_∞ = 1.665e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3030e+02 J
  → Fracture energy : 3.5784e+00 J
  → Total energy    : 2.3388e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=748.16 K
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
  ||ΔD||/||D|| = 2.532e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.948e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=748.16 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.790e-05
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
  ||ΔD||/||D|| = 2.468e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.044e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=748.16 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.395e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.201e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.839e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.385e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=748.16 K
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
  |ΔD|_∞ = 1.664e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3027e+02 J
  → Fracture energy : 3.5874e+00 J
  → Total energy    : 2.3385e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=747.37 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.342e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.213e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.476e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.926e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=747.37 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.711e-05
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
  ||ΔD||/||D|| = 2.417e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.022e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=747.37 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.356e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.152e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.803e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.370e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=747.37 K
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
  ||ΔD||/||D|| = 1.209e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.655e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3024e+02 J
  → Fracture energy : 3.5962e+00 J
  → Total energy    : 2.3383e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=746.60 K
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
  ||ΔD||/||D|| = 2.424e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.883e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=746.60 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.634e-05
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
  ||ΔD||/||D|| = 2.373e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.977e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=746.60 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.317e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.106e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.772e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.336e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=746.60 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.658e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.111e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.188e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.633e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3021e+02 J
  → Fracture energy : 3.6049e+00 J
  → Total energy    : 2.3381e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=745.84 K
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
  |ΔD|_∞ = 2.836e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=745.84 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.558e-05
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
  |ΔD|_∞ = 2.912e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=745.84 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.279e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.061e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.760e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.283e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=745.84 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.640e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.080e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.179e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.597e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3018e+02 J
  → Fracture energy : 3.6136e+00 J
  → Total energy    : 2.3379e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=745.09 K
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
  ||ΔD||/||D|| = 2.404e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.866e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=745.09 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.484e-05
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
  ||ΔD||/||D|| = 2.347e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.951e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=745.09 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.242e-06
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
  ||ΔD||/||D|| = 1.750e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.300e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=745.09 K
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
  |ΔD|_∞ = 1.596e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3014e+02 J
  → Fracture energy : 3.6223e+00 J
  → Total energy    : 2.3377e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=744.34 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.282e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.190e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.371e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.897e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=744.34 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.412e-05
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
  ||ΔD||/||D|| = 2.317e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.985e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=744.34 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.206e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.971e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.729e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.330e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=744.34 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.603e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.019e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.158e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.619e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3011e+02 J
  → Fracture energy : 3.6309e+00 J
  → Total energy    : 2.3375e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=743.60 K
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
  ||ΔD||/||D|| = 2.344e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.912e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=743.60 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.341e-05
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
  ||ΔD||/||D|| = 2.294e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.002e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=743.60 K
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
  |ΔD|_∞ = 2.346e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=743.60 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.585e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.989e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.148e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.632e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3008e+02 J
  → Fracture energy : 3.6394e+00 J
  → Total energy    : 2.3372e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=742.88 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.254e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.177e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.313e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.914e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=742.88 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.272e-05
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
  |ΔD|_∞ = 3.006e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=742.88 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.136e-06
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
  |ΔD|_∞ = 2.351e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=742.88 K
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
  |ΔD|_∞ = 1.638e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3005e+02 J
  → Fracture energy : 3.6479e+00 J
  → Total energy    : 2.3370e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=742.16 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.241e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.169e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.280e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.900e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=742.16 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.205e-05
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
  ||ΔD||/||D|| = 2.233e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.992e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=742.16 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.102e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.836e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.670e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.342e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=742.16 K
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
  |ΔD|_∞ = 1.633e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3002e+02 J
  → Fracture energy : 3.6562e+00 J
  → Total energy    : 2.3368e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=741.45 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.228e-03
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.236e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.863e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=741.45 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.139e-05
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
  |ΔD|_∞ = 2.954e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=741.45 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.069e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.788e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.643e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.314e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=741.45 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.535e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.898e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.104e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.615e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3000e+02 J
  → Fracture energy : 3.6643e+00 J
  → Total energy    : 2.3366e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=740.74 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.215e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.151e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.195e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.798e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=740.74 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.074e-05
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
  |ΔD|_∞ = 2.887e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=740.74 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.037e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.740e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.626e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.263e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=740.74 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.518e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.867e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.093e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.580e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2997e+02 J
  → Fracture energy : 3.6724e+00 J
  → Total energy    : 2.3364e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=740.04 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.202e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.141e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.202e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.736e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=740.04 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.010e-05
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
  |ΔD|_∞ = 2.818e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=740.04 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.005e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.692e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.619e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.195e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=740.04 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.503e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.836e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.087e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.532e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2994e+02 J
  → Fracture energy : 3.6807e+00 J
  → Total energy    : 2.3362e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=739.36 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.190e-03
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.192e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.752e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=739.36 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 5.948e-05
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
  ||ΔD||/||D|| = 2.151e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.837e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=739.36 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.974e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.643e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.611e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.214e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=739.36 K
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
  ||ΔD||/||D|| = 1.083e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.537e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2991e+02 J
  → Fracture energy : 3.6889e+00 J
  → Total energy    : 2.3360e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=738.67 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.177e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.123e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.176e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.760e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=738.67 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 5.887e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.133e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.137e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.849e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=738.67 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.944e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.597e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.601e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.226e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=738.67 K
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
  ||ΔD||/||D|| = 1.077e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.548e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2988e+02 J
  → Fracture energy : 3.6971e+00 J
  → Total energy    : 2.3358e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=738.00 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.166e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.115e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.173e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.755e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=738.00 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 5.828e-05
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
  ||ΔD||/||D|| = 2.130e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.846e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=738.00 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.914e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.552e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.594e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.227e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=738.00 K
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
  ||ΔD||/||D|| = 1.071e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.551e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2985e+02 J
  → Fracture energy : 3.7053e+00 J
  → Total energy    : 2.3355e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=737.33 K
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
  ||ΔD||/||D|| = 2.156e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.739e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=737.33 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 5.769e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.121e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.114e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.832e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=737.33 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.885e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.510e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.583e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.218e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=737.33 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.442e-07
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
  ||ΔD||/||D|| = 1.064e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.546e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2982e+02 J
  → Fracture energy : 3.7135e+00 J
  → Total energy    : 2.3353e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=736.67 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.142e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.102e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.129e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.709e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=736.67 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 5.712e-05
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
  ||ΔD||/||D|| = 2.091e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.802e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=736.67 K
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
  ||ΔD||/||D|| = 1.567e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.197e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=736.67 K
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
  ||ΔD||/||D|| = 1.053e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.533e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2979e+02 J
  → Fracture energy : 3.7215e+00 J
  → Total energy    : 2.3351e+02 J


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
  T_new: min=263.15 K, max=1023.17 K, mean=736.02 K
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
  ||ΔD||/||D|| = 2.105e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.659e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=736.02 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 5.656e-05
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
  ||ΔD||/||D|| = 2.069e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.752e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=736.02 K
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
  ||ΔD||/||D|| = 1.551e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.160e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=736.02 K
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
  ||ΔD||/||D|| = 1.043e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.509e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2976e+02 J
  → Fracture energy : 3.7295e+00 J
  → Total energy    : 2.3349e+02 J


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
  T_new: min=263.15 K, max=1023.18 K, mean=735.37 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.120e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.087e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.096e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.589e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=735.37 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 5.601e-05
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
  |ΔD|_∞ = 2.681e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=735.37 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.800e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.392e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.542e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.106e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=735.37 K
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
  |ΔD|_∞ = 1.473e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2973e+02 J
  → Fracture energy : 3.7376e+00 J
  → Total energy    : 2.3347e+02 J

Simulation completed in 1722.26 s
Total time steps solved: 100
