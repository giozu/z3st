Info    : Reading 'mesh.msh'...
Info    : 34152 nodes
Info    : 68322 elements
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
    Set 34152 DOFs to 1023.15 K
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
  T_new: min=263.15 K, max=1043.42 K, mean=993.46 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.281e-01
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
  T_new: min=263.15 K, max=1043.42 K, mean=993.46 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.009e-03
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
  ||ΔD||/||D|| = 6.295e-01
  [adaptive] relax_D=0.55
  |ΔD|_∞ = 3.080e-01

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1043.42 K, mean=993.46 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.542e-03
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
  ||ΔD||/||D|| = 1.599e+00
  [adaptive] relax_D=0.52
  |ΔD|_∞ = 4.243e-01

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1043.42 K, mean=993.46 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.997e-04
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.457e+00
  [adaptive] relax_D=0.50
  |ΔD|_∞ = 2.116e-01

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1043.42 K, mean=993.46 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.987e-06
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
  ||ΔD||/||D|| = 1.989e+00
  [adaptive] relax_D=0.47
  |ΔD|_∞ = 1.088e-01

Convergence check


#### Iteration 6/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1043.42 K, mean=993.46 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.994e-07
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
  ||ΔD||/||D|| = 1.145e+00
  [adaptive] relax_D=0.52
  |ΔD|_∞ = 5.802e-02

Convergence check


#### Iteration 7/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1043.42 K, mean=993.46 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.497e-08
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
  ||ΔD||/||D|| = 6.938e-01
  [adaptive] relax_D=0.57
  |ΔD|_∞ = 3.457e-02

Convergence check


#### Iteration 8/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1043.42 K, mean=993.46 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.248e-09
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
  ||ΔD||/||D|| = 3.689e-01
  [adaptive] relax_D=0.63
  |ΔD|_∞ = 1.835e-02

Convergence check


#### Iteration 9/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1043.42 K, mean=993.46 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.242e-11
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
  ||ΔD||/||D|| = 1.743e-01
  [adaptive] relax_D=0.69
  |ΔD|_∞ = 8.671e-03

Convergence check


#### Iteration 10/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1043.42 K, mean=993.46 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.121e-12
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.141e-02
  [adaptive] relax_D=0.76
  |ΔD|_∞ = 3.552e-03

Convergence check


#### Iteration 11/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1043.42 K, mean=993.46 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.560e-13
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
  ||ΔD||/||D|| = 2.432e-02
  [adaptive] relax_D=0.84
  |ΔD|_∞ = 1.210e-03

Convergence check


#### Iteration 12/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1043.42 K, mean=993.46 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.803e-15
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
  ||ΔD||/||D|| = 6.435e-03
  [adaptive] relax_D=0.92
  |ΔD|_∞ = 3.200e-04

Convergence check


#### Iteration 13/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1043.42 K, mean=993.46 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.901e-16
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.026e-12
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.165e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.795e-05

Convergence check


#### Iteration 14/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1043.42 K, mean=993.46 K
  T^n (self.T): min=1023.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.280e-17
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.497e-13
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.765e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.857e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 14 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3478e+02 J
  → Fracture energy : 7.9499e-01 J
  → Total energy    : 2.3557e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=953.39 K
  T^n (self.T): min=263.15 K, max=1043.42 K
  ||ΔT||/||T|| = 9.708e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.363e-03
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.348e-01
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.878e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=953.39 K
  T^n (self.T): min=263.15 K, max=1043.42 K
  ||ΔT||/||T|| = 4.854e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.469e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.995e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.644e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=953.39 K
  T^n (self.T): min=263.15 K, max=1043.42 K
  ||ΔT||/||T|| = 2.427e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.134e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.190e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.304e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=953.39 K
  T^n (self.T): min=263.15 K, max=1043.42 K
  ||ΔT||/||T|| = 1.214e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.672e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.888e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.645e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=953.39 K
  T^n (self.T): min=263.15 K, max=1043.42 K
  ||ΔT||/||T|| = 6.068e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.838e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.092e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.967e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3247e+02 J
  → Fracture energy : 2.3959e+00 J
  → Total energy    : 2.3487e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=930.43 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.599e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.497e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.497e-01
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.148e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=930.43 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.299e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.885e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.894e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.984e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=930.43 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.150e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.981e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.115e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.911e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=930.43 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.749e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.574e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.262e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.565e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=930.43 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.874e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.338e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.651e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.162e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3261e+02 J
  → Fracture energy : 3.2653e+00 J
  → Total energy    : 2.3587e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=913.79 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.094e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.450e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.493e-01
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.679e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=913.79 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.547e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.487e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.095e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.587e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=913.79 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.734e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.499e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.177e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.492e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=913.79 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.867e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.138e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.487e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.206e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=913.79 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.934e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.015e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.116e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.875e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3234e+02 J
  → Fracture energy : 3.8852e+00 J
  → Total energy    : 2.3622e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=900.50 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.370e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.317e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.077e-01
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.455e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=900.50 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.185e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.130e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.771e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.423e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=900.50 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.926e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.110e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.324e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.293e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=900.50 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.963e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.833e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.491e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.020e-04

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=900.50 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.481e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.807e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.550e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.386e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3202e+02 J
  → Fracture energy : 4.3945e+00 J
  → Total energy    : 2.3641e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=889.32 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.939e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.171e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.428e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.355e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=889.32 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.696e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.809e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.048e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.339e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=889.32 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.848e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.810e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.388e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.178e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=889.32 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.424e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.612e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.004e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.155e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=889.32 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.212e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.660e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.269e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.555e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3171e+02 J
  → Fracture energy : 4.8374e+00 J
  → Total energy    : 2.3655e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=879.63 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.650e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.106e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.906e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.352e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=879.63 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 8.250e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.575e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.971e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.277e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=879.63 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.125e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.590e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.863e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.109e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=879.63 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.063e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.450e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.749e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.480e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=879.63 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.031e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.554e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.119e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.009e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3142e+02 J
  → Fracture energy : 5.2325e+00 J
  → Total energy    : 2.3665e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=871.04 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.441e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.129e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.889e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.421e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=871.04 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.207e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.412e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.352e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.222e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=871.04 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.604e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.427e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.660e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.050e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=871.04 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.802e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.329e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.666e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.980e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=871.04 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.009e-08
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.474e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.061e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.618e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3115e+02 J
  → Fracture energy : 5.5926e+00 J
  → Total energy    : 2.3674e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=863.30 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.283e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.229e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.349e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.518e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=863.30 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.415e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.294e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.203e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.271e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=863.30 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.208e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.299e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.713e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.928e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=863.30 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.604e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.233e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.705e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.489e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=863.30 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 8.019e-08
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.409e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.068e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.256e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3087e+02 J
  → Fracture energy : 5.9307e+00 J
  → Total energy    : 2.3680e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=856.26 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.158e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.459e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.195e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.614e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=856.26 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.792e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.210e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.295e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.392e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=856.26 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.896e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.188e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.864e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.509e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=856.26 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.448e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.147e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.802e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.039e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=856.26 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.240e-08
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.353e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.111e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.929e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3056e+02 J
  → Fracture energy : 6.2621e+00 J
  → Total energy    : 2.3682e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=849.78 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.057e-02
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.823e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.073e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.650e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=849.78 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.287e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.138e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.306e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.448e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=849.78 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.643e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.074e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.927e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.014e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=849.78 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.322e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.059e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.850e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.527e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=849.78 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.609e-08
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.295e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.134e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.568e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.3024e+02 J
  → Fracture energy : 6.5868e+00 J
  → Total energy    : 2.3682e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=843.79 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.737e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.628e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.801e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.624e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=843.79 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.869e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.181e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.118e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.438e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=843.79 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.434e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.986e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.828e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.012e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=843.79 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.217e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.976e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.794e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.483e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=843.79 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 6.086e-08
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.237e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.099e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.190e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2991e+02 J
  → Fracture energy : 6.8913e+00 J
  → Total energy    : 2.3680e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=838.20 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.032e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.712e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.369e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.529e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=838.20 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.516e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.356e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.758e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.374e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=838.20 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.258e-05
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
  ||ΔD||/||D|| = 2.590e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.826e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=838.20 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.129e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.907e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.647e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.406e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=838.20 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.645e-08
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.184e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.009e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.950e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2959e+02 J
  → Fracture energy : 7.1645e+00 J
  → Total energy    : 2.3675e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=832.96 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 8.429e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.342e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.839e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.441e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=832.96 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 4.214e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.063e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.304e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.293e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=832.96 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 2.107e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.749e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.280e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.205e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=832.96 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.054e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.789e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.452e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.060e-05

Convergence check


#### Iteration 5/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=832.96 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 5.268e-08
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.115e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.907e-07
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.771e-06

Convergence check

**[SUCCESS]** Staggered solver converged in 5 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2930e+02 J
  → Fracture energy : 7.4019e+00 J
  → Total energy    : 2.3670e+02 J
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
  T_new: min=263.15 K, max=1023.15 K, mean=828.04 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.906e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 5.493e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.317e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.443e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=828.04 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.953e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.652e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.860e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.289e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=828.04 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.977e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.540e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.974e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.212e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.15 K, mean=828.04 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.883e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.676e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.258e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.972e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2904e+02 J
  → Fracture energy : 7.6060e+00 J
  → Total energy    : 2.3665e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=823.38 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 7.449e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 4.109e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.870e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.459e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=823.38 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 3.724e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.133e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.478e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.318e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=823.38 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 1.862e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.322e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.711e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.425e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=823.38 K
  T^n (self.T): min=263.15 K, max=1023.15 K
  ||ΔT||/||T|| = 9.311e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.571e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.092e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.164e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2881e+02 J
  → Fracture energy : 7.7829e+00 J
  → Total energy    : 2.3659e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=818.97 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.045e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 3.218e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.526e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.458e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=818.97 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.522e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.840e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.189e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.328e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=818.97 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.761e-05
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
  ||ΔD||/||D|| = 1.511e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.553e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=818.97 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.806e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.503e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.633e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.228e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2859e+02 J
  → Fracture energy : 7.9390e+00 J
  → Total energy    : 2.3653e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=814.78 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.685e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.975e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.286e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.499e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=814.78 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.342e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.744e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.987e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.327e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=814.78 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.671e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.131e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.372e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.619e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=814.78 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.356e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.460e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.731e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.292e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2839e+02 J
  → Fracture energy : 8.0795e+00 J
  → Total energy    : 2.3647e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=810.78 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.362e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.905e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.113e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.577e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=810.78 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.181e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.689e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.841e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.402e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=810.78 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.591e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.082e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.272e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.873e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=810.78 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.953e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.423e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.088e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.313e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2820e+02 J
  → Fracture energy : 8.2085e+00 J
  → Total energy    : 2.3641e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=806.96 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.071e-03
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
  ||ΔD||/||D|| = 1.960e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.622e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=806.96 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.036e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.647e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.715e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.446e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=806.96 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.518e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.037e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.187e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.027e-03

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=806.96 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.589e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.389e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.541e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.607e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2801e+02 J
  → Fracture energy : 8.3278e+00 J
  → Total energy    : 2.3634e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=803.31 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.807e-03
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.813e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.566e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=803.31 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.904e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.597e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.594e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.410e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=803.31 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.452e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.992e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.104e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.992e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=803.31 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.259e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.356e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.001e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.464e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2783e+02 J
  → Fracture energy : 8.4388e+00 J
  → Total energy    : 2.3627e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=799.80 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.567e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.776e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.685e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.503e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=799.80 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.783e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.553e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.487e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.363e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=799.80 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.392e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.954e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.032e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.745e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=799.80 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.958e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.328e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.544e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.344e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2766e+02 J
  → Fracture energy : 8.5429e+00 J
  → Total energy    : 2.3620e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=796.43 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.346e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.721e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.564e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.411e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=796.43 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.673e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.513e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.389e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.287e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=796.43 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.337e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.920e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.684e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.278e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=796.43 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.683e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.302e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.163e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.088e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2749e+02 J
  → Fracture energy : 8.6409e+00 J
  → Total energy    : 2.3613e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=793.19 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.144e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.684e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.468e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.316e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=793.19 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.572e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.481e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.310e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.209e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=793.19 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.286e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.890e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.170e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.755e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=793.19 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.430e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.280e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.853e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.734e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2732e+02 J
  → Fracture energy : 8.7345e+00 J
  → Total energy    : 2.3606e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=790.07 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.957e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.671e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.391e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.225e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=790.07 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.478e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.455e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.251e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.133e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=790.07 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.239e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.863e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.800e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.295e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=790.07 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.196e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.259e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.630e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.489e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2716e+02 J
  → Fracture energy : 8.8244e+00 J
  → Total energy    : 2.3598e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=787.05 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.784e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.660e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.339e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.142e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=787.05 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.392e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.432e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.208e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.062e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=787.05 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.196e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.839e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.518e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.783e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=787.05 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.980e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.240e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.464e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.150e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2700e+02 J
  → Fracture energy : 8.9109e+00 J
  → Total energy    : 2.3591e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=784.13 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.623e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.662e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.292e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.075e-01

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=784.13 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.311e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.412e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.175e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.011e-02

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=784.13 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.156e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.816e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.327e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.467e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=784.13 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.779e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.222e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.354e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.966e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2684e+02 J
  → Fracture energy : 8.9951e+00 J
  → Total energy    : 2.3584e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=781.31 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.473e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.676e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.262e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.779e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=781.31 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.237e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.395e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.152e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.303e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=781.31 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.118e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.795e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.197e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.878e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=781.31 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.592e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.205e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.286e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.569e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2668e+02 J
  → Fracture energy : 9.0775e+00 J
  → Total energy    : 2.3576e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=778.58 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.334e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.694e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.230e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.171e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=778.58 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.167e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.379e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.133e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.812e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=778.58 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.083e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.775e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.094e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.661e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=778.58 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.417e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.189e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.234e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.519e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2653e+02 J
  → Fracture energy : 9.1589e+00 J
  → Total energy    : 2.3569e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=775.93 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.203e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.696e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.202e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 9.114e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=775.93 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.101e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.356e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.111e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.772e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=775.93 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.051e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.752e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.971e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.588e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=775.93 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.253e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.172e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.168e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.452e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2637e+02 J
  → Fracture energy : 9.2388e+00 J
  → Total energy    : 2.3561e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=773.36 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.080e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.681e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.168e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.947e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=773.36 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.040e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.331e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.086e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.654e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=773.36 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.020e-05
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.819e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.482e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=773.36 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 5.100e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.156e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.085e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.391e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2622e+02 J
  → Fracture energy : 9.3170e+00 J
  → Total energy    : 2.3554e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=770.86 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.965e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.650e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.135e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.806e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=770.86 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.982e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.302e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.059e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.529e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=770.86 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 9.911e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.706e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.648e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.434e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=770.86 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.956e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.140e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.988e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.361e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2607e+02 J
  → Fracture energy : 9.3933e+00 J
  → Total energy    : 2.3546e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=768.43 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.856e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.615e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.097e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.381e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=768.43 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.928e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.273e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.029e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.185e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=768.43 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 9.639e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.683e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.467e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.176e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=768.43 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.820e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.124e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.884e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.177e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2592e+02 J
  → Fracture energy : 9.4680e+00 J
  → Total energy    : 2.3539e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=766.07 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.753e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.582e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.060e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.332e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=766.07 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.877e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.243e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.986e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.080e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=766.07 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 9.383e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.660e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.269e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.072e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=766.07 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.691e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.108e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.767e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.095e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2577e+02 J
  → Fracture energy : 9.5408e+00 J
  → Total energy    : 2.3532e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=763.77 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.656e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.545e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.021e-02
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.227e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=763.77 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.828e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.213e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.647e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.960e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=763.77 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 9.140e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.638e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.042e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.006e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=763.77 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.570e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.093e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.631e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.074e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2563e+02 J
  → Fracture energy : 9.6114e+00 J
  → Total energy    : 2.3524e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=761.53 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.564e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.511e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.830e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 8.079e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=761.53 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.782e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.185e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.335e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.888e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=761.53 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.911e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.617e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.833e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.955e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=761.53 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.455e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.079e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.502e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.030e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2549e+02 J
  → Fracture energy : 9.6807e+00 J
  → Total energy    : 2.3517e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=759.35 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.477e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.472e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.442e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.946e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=759.35 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.739e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.156e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.986e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.748e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=759.35 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.693e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.596e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.593e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.875e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=759.35 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.346e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.065e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.353e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.000e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2535e+02 J
  → Fracture energy : 9.7485e+00 J
  → Total energy    : 2.3509e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=757.22 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.394e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.418e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 9.022e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.600e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=757.22 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.697e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.124e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.622e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.405e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=757.22 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.486e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.574e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.344e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.611e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=757.22 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.243e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.051e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.198e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.806e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2521e+02 J
  → Fracture energy : 9.8142e+00 J
  → Total energy    : 2.3502e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=755.14 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.316e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.361e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.647e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.521e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=755.14 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.658e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.091e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.272e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.337e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=755.14 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.289e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.553e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.096e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.589e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=755.14 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.144e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.037e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.040e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.816e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2507e+02 J
  → Fracture energy : 9.8780e+00 J
  → Total energy    : 2.3495e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=753.10 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.240e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.311e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 8.263e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.311e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=753.10 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.620e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.061e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.954e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.098e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=753.10 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 8.101e-06
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
  ||ΔD||/||D|| = 5.880e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.388e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=753.10 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 4.051e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.024e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.906e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.695e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2494e+02 J
  → Fracture energy : 9.9399e+00 J
  → Total energy    : 2.3488e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=751.12 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.169e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.266e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.984e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.293e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=751.12 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.584e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.033e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.683e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.115e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=751.12 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.922e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.513e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.686e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.352e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=751.12 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.961e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.011e-07
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.782e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.607e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2480e+02 J
  → Fracture energy : 1.0000e+01 J
  → Total energy    : 2.3480e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=749.18 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.101e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.225e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.687e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.089e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=749.18 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.550e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.006e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.417e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.925e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=749.18 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.752e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.495e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.498e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.237e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=749.18 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.876e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.994e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.663e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.552e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2467e+02 J
  → Fracture energy : 1.0059e+01 J
  → Total energy    : 2.3473e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=747.28 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.035e-03
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
  ||ΔD||/||D|| = 7.417e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.093e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=747.28 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.518e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.981e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.182e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.939e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=747.28 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.588e-06
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
  ||ΔD||/||D|| = 5.338e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.216e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=747.28 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.794e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.880e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.563e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.509e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2454e+02 J
  → Fracture energy : 1.0116e+01 J
  → Total energy    : 2.3466e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=745.42 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.973e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.150e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 7.167e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.107e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=745.42 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.486e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.958e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.942e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.970e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=745.42 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.432e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.461e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.169e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.270e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=745.42 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.716e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.773e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.456e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.570e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2442e+02 J
  → Fracture energy : 1.0172e+01 J
  → Total energy    : 2.3459e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=743.60 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.913e-03
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.917e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 7.025e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=743.60 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.457e-04
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
  ||ΔD||/||D|| = 6.731e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.885e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=743.60 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.283e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.446e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.023e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.225e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=743.60 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.641e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.672e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.364e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.556e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2429e+02 J
  → Fracture energy : 1.0226e+01 J
  → Total energy    : 2.3452e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=741.82 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.856e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.099e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.720e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.746e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=741.82 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.428e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.918e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.548e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.578e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=741.82 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.139e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.432e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.892e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.996e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=741.82 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.570e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.574e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.280e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.409e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2417e+02 J
  → Fracture energy : 1.0279e+01 J
  → Total energy    : 2.3445e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=740.07 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.801e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.070e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.519e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.458e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=740.07 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.400e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.898e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.355e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.359e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=740.07 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 7.002e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.417e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.759e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.831e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=740.07 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.501e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.476e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.198e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.288e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2405e+02 J
  → Fracture energy : 1.0330e+01 J
  → Total energy    : 2.3438e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=738.36 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.748e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.039e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.300e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.435e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=738.36 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.374e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.877e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.175e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.361e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=738.36 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.869e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.403e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.636e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.818e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=738.36 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.435e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.382e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.121e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.264e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2393e+02 J
  → Fracture energy : 1.0379e+01 J
  → Total energy    : 2.3430e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=736.68 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.697e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 2.011e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.136e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.442e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=736.68 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.348e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.858e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 6.013e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.390e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=736.68 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.742e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.389e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.516e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.867e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=736.68 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.371e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.292e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.043e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.316e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2381e+02 J
  → Fracture energy : 1.0428e+01 J
  → Total energy    : 2.3423e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=735.03 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.648e-03
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
  ||ΔD||/||D|| = 5.934e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.262e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=735.03 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.324e-04
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
  ||ΔD||/||D|| = 5.831e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 6.224e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=735.03 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.620e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.376e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.388e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.760e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=735.03 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.310e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.205e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.961e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.259e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2369e+02 J
  → Fracture energy : 1.0476e+01 J
  → Total energy    : 2.3416e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=733.42 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.601e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.946e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.722e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.938e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=733.42 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.300e-04
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
  ||ΔD||/||D|| = 5.639e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.900e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=733.42 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.502e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.363e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.252e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.525e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=733.42 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.251e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.120e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.875e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.108e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2357e+02 J
  → Fracture energy : 1.0522e+01 J
  → Total energy    : 2.3409e+02 J
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
  T_new: min=263.15 K, max=1023.16 K, mean=731.83 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.556e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.913e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.540e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.652e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=731.83 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.278e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.800e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.463e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.650e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=731.83 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.389e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.350e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.126e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.329e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.16 K, mean=731.83 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.194e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 9.036e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.793e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.965e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2346e+02 J
  → Fracture energy : 1.0567e+01 J
  → Total energy    : 2.3402e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=730.27 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 2.512e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.883e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.362e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.477e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=730.27 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 1.256e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.782e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.295e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.481e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=730.27 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 6.280e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.337e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.004e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.213e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=730.27 K
  T^n (self.T): min=263.15 K, max=1023.16 K
  ||ΔT||/||T|| = 3.140e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.954e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.714e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.896e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2334e+02 J
  → Fracture energy : 1.0612e+01 J
  → Total energy    : 2.3396e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=728.74 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.470e-03
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
  ||ΔD||/||D|| = 5.199e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.241e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=728.74 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.235e-04
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
  ||ΔD||/||D|| = 5.147e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.245e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=728.74 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.174e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.325e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.898e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.041e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=728.74 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.087e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.875e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.644e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.786e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2323e+02 J
  → Fracture energy : 1.0655e+01 J
  → Total energy    : 2.3389e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=727.23 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.429e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.833e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.049e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.067e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=727.23 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.214e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.748e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 5.005e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 5.094e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=727.23 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 6.072e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.314e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.795e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.912e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=727.23 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.036e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.798e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.578e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.709e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2312e+02 J
  → Fracture energy : 1.0697e+01 J
  → Total energy    : 2.3382e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=725.76 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.389e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.807e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.895e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.931e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=725.76 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.195e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.731e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.861e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.972e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=725.76 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 5.974e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.302e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.689e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.843e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=725.76 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.987e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.725e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.508e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.667e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2301e+02 J
  → Fracture energy : 1.0739e+01 J
  → Total energy    : 2.3375e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=724.30 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.351e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.780e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.759e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.737e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=724.30 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.176e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.715e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.728e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.784e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=724.30 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 5.878e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.291e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.591e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.715e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=724.30 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.939e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.653e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.443e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.579e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2290e+02 J
  → Fracture energy : 1.0780e+01 J
  → Total energy    : 2.3368e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=722.87 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.315e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.756e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.627e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.617e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=722.87 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.157e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.700e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.602e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.718e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=722.87 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 5.786e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.281e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.498e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.678e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=722.87 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.893e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.584e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.381e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.557e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2279e+02 J
  → Fracture energy : 1.0820e+01 J
  → Total energy    : 2.3361e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=721.47 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.279e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.736e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.502e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.587e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=721.47 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.139e-04
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
  ||ΔD||/||D|| = 4.478e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.685e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=721.47 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 5.697e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.271e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.406e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.655e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=721.47 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.849e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.516e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.320e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.544e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2269e+02 J
  → Fracture energy : 1.0859e+01 J
  → Total energy    : 2.3355e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=720.08 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.244e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.718e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.378e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.487e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=720.08 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.122e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.672e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.366e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.579e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=720.08 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 5.611e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.261e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.324e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.575e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=720.08 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.805e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.451e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.266e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.492e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2258e+02 J
  → Fracture energy : 1.0898e+01 J
  → Total energy    : 2.3348e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=718.72 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.211e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.699e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.265e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.321e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=718.72 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.105e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.658e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.250e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.406e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=718.72 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 5.527e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.251e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.236e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.442e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=718.72 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.764e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.387e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.207e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.401e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2248e+02 J
  → Fracture energy : 1.0936e+01 J
  → Total energy    : 2.3341e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=717.38 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.178e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.679e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.143e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.095e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=717.38 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.089e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.644e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.134e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.222e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=717.38 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 5.446e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.241e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.150e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.310e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=717.38 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.723e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.324e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.150e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.312e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2237e+02 J
  → Fracture energy : 1.0973e+01 J
  → Total energy    : 2.3334e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=716.06 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.147e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.659e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.023e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.046e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=716.06 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.073e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.631e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 4.017e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.176e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=716.06 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 5.367e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.232e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.062e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.279e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=716.06 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.684e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.264e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.090e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.294e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2227e+02 J
  → Fracture energy : 1.1010e+01 J
  → Total energy    : 2.3328e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=714.76 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.116e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.640e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.905e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.966e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=714.76 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.058e-04
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.900e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 4.095e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=714.76 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 5.291e-06
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.974e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.220e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=714.76 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.646e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.205e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.030e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.256e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2217e+02 J
  → Fracture energy : 1.1046e+01 J
  → Total energy    : 2.3321e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=713.48 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.087e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.622e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.800e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.838e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=713.48 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.043e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.606e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.798e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.965e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=713.48 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 5.217e-06
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
  ||ΔD||/||D|| = 2.897e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.122e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=713.48 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.608e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.147e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.979e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.191e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2207e+02 J
  → Fracture energy : 1.1081e+01 J
  → Total energy    : 2.3315e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=712.22 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.058e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.605e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.732e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.665e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=712.22 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.029e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.593e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.724e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.789e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=712.22 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 5.145e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.206e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.839e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.987e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=712.22 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.572e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.090e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.939e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.099e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2197e+02 J
  → Fracture energy : 1.1116e+01 J
  → Total energy    : 2.3308e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=710.97 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.030e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.588e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.659e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.511e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=710.97 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.015e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.581e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.652e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.618e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=710.97 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 5.075e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.197e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.783e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.854e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=710.97 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.537e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 8.032e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.900e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.007e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2187e+02 J
  → Fracture energy : 1.1150e+01 J
  → Total energy    : 2.3302e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=709.75 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.003e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.572e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.578e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.446e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=709.75 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.001e-04
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.569e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.573e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.587e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=709.75 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 5.007e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.188e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.725e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.841e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=709.75 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.503e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.977e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.861e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.002e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2177e+02 J
  → Fracture energy : 1.1184e+01 J
  → Total energy    : 2.3295e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=708.54 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.976e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.557e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.493e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.393e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=708.54 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 9.881e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.558e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.489e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.548e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=708.54 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.941e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.180e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.661e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.813e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=708.54 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.470e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.923e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.818e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.984e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2167e+02 J
  → Fracture energy : 1.1217e+01 J
  → Total energy    : 2.3289e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=707.34 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.951e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.542e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.406e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.329e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=707.34 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 9.753e-05
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
  ||ΔD||/||D|| = 3.402e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.483e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=707.34 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.876e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.172e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.596e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.764e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=707.34 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.438e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.870e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.775e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.952e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2157e+02 J
  → Fracture energy : 1.1250e+01 J
  → Total energy    : 2.3282e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=706.17 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.925e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.527e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.318e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.271e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=706.17 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 9.627e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.535e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.316e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.388e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=706.17 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.814e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.164e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.532e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.691e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=706.17 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.407e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.818e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.733e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.902e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2148e+02 J
  → Fracture energy : 1.1282e+01 J
  → Total energy    : 2.3276e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=705.01 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.901e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.511e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.251e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.199e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=705.01 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 9.505e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.524e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.255e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.311e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=705.01 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.753e-06
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
  ||ΔD||/||D|| = 2.487e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.605e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=705.01 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.376e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.768e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.702e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.836e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2138e+02 J
  → Fracture energy : 1.1314e+01 J
  → Total energy    : 2.3270e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=703.87 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.877e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.496e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.209e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.101e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=703.87 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 9.386e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.513e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.206e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.209e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=703.87 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.693e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.149e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.448e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.528e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=703.87 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.347e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.717e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.674e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.773e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2129e+02 J
  → Fracture energy : 1.1346e+01 J
  → Total energy    : 2.3263e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=702.74 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.854e-03
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.156e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.072e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=702.74 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 9.271e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.503e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.153e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.165e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=702.74 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.635e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.141e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.407e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.480e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=702.74 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.318e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.668e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.646e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.757e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2119e+02 J
  → Fracture energy : 1.1377e+01 J
  → Total energy    : 2.3257e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=701.62 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.832e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.472e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.106e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.107e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=701.62 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 9.158e-05
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
  ||ΔD||/||D|| = 3.103e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.202e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=701.62 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.579e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.134e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.369e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.501e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=701.62 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.289e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.621e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.620e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.740e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2110e+02 J
  → Fracture energy : 1.1408e+01 J
  → Total energy    : 2.3251e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=700.52 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.810e-03
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
  ||ΔD||/||D|| = 3.059e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.112e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=700.52 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 9.048e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.484e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.052e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.207e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=700.52 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.524e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.127e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.329e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.508e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=700.52 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.262e-07
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
  ||ΔD||/||D|| = 1.592e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.747e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2101e+02 J
  → Fracture energy : 1.1438e+01 J
  → Total energy    : 2.3245e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=699.44 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.788e-03
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 3.000e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.081e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=699.44 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 8.940e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.475e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.996e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.176e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=699.44 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.470e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.121e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.287e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.486e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=699.44 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.235e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.531e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.565e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.734e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2092e+02 J
  → Fracture energy : 1.1468e+01 J
  → Total energy    : 2.3238e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=698.36 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.767e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.445e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.941e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.015e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=698.36 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 8.835e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.467e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.941e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.109e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=698.36 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.418e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.115e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.246e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.437e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=698.36 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.209e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.489e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.536e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.701e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2082e+02 J
  → Fracture energy : 1.1498e+01 J
  → Total energy    : 2.3232e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=697.31 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.747e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.435e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.896e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.986e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=697.31 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 8.733e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.459e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.892e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.084e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=697.31 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.366e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.108e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.207e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.410e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=697.31 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.183e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.448e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.509e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.676e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2073e+02 J
  → Fracture energy : 1.1527e+01 J
  → Total energy    : 2.3226e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=696.26 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.727e-03
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
  ||ΔD||/||D|| = 2.851e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.988e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=696.26 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 8.633e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.450e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.847e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.091e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=696.26 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.317e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.102e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.172e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.421e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=696.26 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.158e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.407e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.486e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.687e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2064e+02 J
  → Fracture energy : 1.1557e+01 J
  → Total energy    : 2.3220e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=695.23 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.707e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.414e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.805e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.959e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=695.23 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 8.536e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.442e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.798e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.066e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=695.23 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.268e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.096e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.134e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.406e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=695.23 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.134e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.367e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.459e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.680e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2055e+02 J
  → Fracture energy : 1.1585e+01 J
  → Total energy    : 2.3214e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=694.21 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.688e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.404e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.751e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.916e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=694.21 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 8.440e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.434e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.745e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 3.027e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=694.21 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.220e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.090e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.093e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.379e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=694.21 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.110e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.328e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.431e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.664e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2046e+02 J
  → Fracture energy : 1.1614e+01 J
  → Total energy    : 2.3208e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=693.20 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.669e-03
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
  ||ΔD||/||D|| = 2.697e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.875e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=693.20 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 8.347e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.426e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.691e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.986e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=693.20 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.174e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.085e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.052e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.351e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=693.20 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.087e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.289e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.403e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.647e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2038e+02 J
  → Fracture energy : 1.1642e+01 J
  → Total energy    : 2.3202e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=692.20 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.651e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.385e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.652e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.824e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=692.20 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 8.256e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.418e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.647e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.934e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=692.20 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.128e-06
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
  ||ΔD||/||D|| = 2.018e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.312e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=692.20 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.064e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.251e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.379e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.622e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2029e+02 J
  → Fracture energy : 1.1669e+01 J
  → Total energy    : 2.3196e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=691.22 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.633e-03
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
  ||ΔD||/||D|| = 2.618e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.762e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=691.22 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 8.167e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.411e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.607e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.870e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=691.22 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.084e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.073e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.986e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.263e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=691.22 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.042e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.213e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.357e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.589e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2020e+02 J
  → Fracture energy : 1.1697e+01 J
  → Total energy    : 2.3190e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=690.24 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.616e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.367e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.571e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.686e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=690.24 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 8.080e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.403e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.564e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.790e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=690.24 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 4.040e-06
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
  ||ΔD||/||D|| = 1.954e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.202e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=690.24 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 2.020e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.176e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.336e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.547e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2012e+02 J
  → Fracture energy : 1.1723e+01 J
  → Total energy    : 2.3184e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=689.28 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.599e-03
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
  ||ΔD||/||D|| = 2.531e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.596e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=689.28 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 7.995e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.396e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.525e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.695e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=689.28 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.997e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.062e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.926e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.127e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=689.28 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.999e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.141e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.317e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.495e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2003e+02 J
  → Fracture energy : 1.1750e+01 J
  → Total energy    : 2.3178e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=688.33 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.582e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.354e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.489e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.546e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=688.33 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 7.911e-05
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.489e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.644e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=688.33 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.956e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.057e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.900e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.105e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=688.33 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.978e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.107e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.300e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.489e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1994e+02 J
  → Fracture energy : 1.1776e+01 J
  → Total energy    : 2.3172e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=687.39 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.566e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.348e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.462e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.543e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=687.39 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 7.830e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.383e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.459e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.680e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=687.39 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.915e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.052e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.876e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.133e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=687.39 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.957e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.073e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.284e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.508e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1986e+02 J
  → Fracture energy : 1.1802e+01 J
  → Total energy    : 2.3166e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=686.46 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.550e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.342e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.439e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.580e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=686.46 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 7.750e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.377e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.437e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.715e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=686.46 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.875e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.048e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.860e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.159e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=686.46 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.938e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.040e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.273e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.525e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1978e+02 J
  → Fracture energy : 1.1828e+01 J
  → Total energy    : 2.3160e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=685.54 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.534e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.337e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.413e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.614e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=685.54 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 7.672e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.371e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.416e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.745e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=685.54 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.836e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.043e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.845e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.181e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=685.54 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.918e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 7.007e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.264e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.540e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1969e+02 J
  → Fracture energy : 1.1854e+01 J
  → Total energy    : 2.3155e+02 J
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
  T_new: min=263.15 K, max=1023.17 K, mean=684.63 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.519e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.331e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.401e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.645e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=684.63 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 7.595e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.365e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.400e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.770e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=684.63 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.798e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.038e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.832e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.198e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.17 K, mean=684.63 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.899e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.975e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.254e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.550e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1961e+02 J
  → Fracture energy : 1.1880e+01 J
  → Total energy    : 2.3149e+02 J
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
  T_new: min=263.15 K, max=1023.18 K, mean=683.73 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.504e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.326e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.379e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.668e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=683.73 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 7.520e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.358e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.385e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.786e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=683.73 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 3.760e-06
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.823e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.207e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=683.73 K
  T^n (self.T): min=263.15 K, max=1023.17 K
  ||ΔT||/||T|| = 1.880e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.943e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.248e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.555e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1952e+02 J
  → Fracture energy : 1.1906e+01 J
  → Total energy    : 2.3143e+02 J
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
  T_new: min=263.15 K, max=1023.18 K, mean=682.84 K
  T^n (self.T): min=263.15 K, max=1023.18 K
  ||ΔT||/||T|| = 1.489e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.321e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.372e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.663e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=682.84 K
  T^n (self.T): min=263.15 K, max=1023.18 K
  ||ΔT||/||T|| = 7.447e-05
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
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.372e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.774e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=682.84 K
  T^n (self.T): min=263.15 K, max=1023.18 K
  ||ΔT||/||T|| = 3.723e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.028e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.812e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.194e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=682.84 K
  T^n (self.T): min=263.15 K, max=1023.18 K
  ||ΔT||/||T|| = 1.862e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.911e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.240e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.544e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1944e+02 J
  → Fracture energy : 1.1931e+01 J
  → Total energy    : 2.3137e+02 J
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
  T_new: min=263.15 K, max=1023.18 K, mean=681.95 K
  T^n (self.T): min=263.15 K, max=1023.18 K
  ||ΔT||/||T|| = 1.475e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.316e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.355e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.619e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=681.95 K
  T^n (self.T): min=263.15 K, max=1023.18 K
  ||ΔT||/||T|| = 7.375e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.347e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.359e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.722e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=681.95 K
  T^n (self.T): min=263.15 K, max=1023.18 K
  ||ΔT||/||T|| = 3.687e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.024e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.804e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.150e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=681.95 K
  T^n (self.T): min=263.15 K, max=1023.18 K
  ||ΔT||/||T|| = 1.844e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.880e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.235e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.513e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1936e+02 J
  → Fracture energy : 1.1957e+01 J
  → Total energy    : 2.3131e+02 J
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
  T_new: min=263.15 K, max=1023.18 K, mean=681.08 K
  T^n (self.T): min=263.15 K, max=1023.18 K
  ||ΔT||/||T|| = 1.461e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.313e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.333e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.577e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=681.08 K
  T^n (self.T): min=263.15 K, max=1023.18 K
  ||ΔT||/||T|| = 7.304e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.341e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.338e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.703e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=681.08 K
  T^n (self.T): min=263.15 K, max=1023.18 K
  ||ΔT||/||T|| = 3.652e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.020e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.789e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.141e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=681.08 K
  T^n (self.T): min=263.15 K, max=1023.18 K
  ||ΔT||/||T|| = 1.826e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.850e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.226e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.508e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1928e+02 J
  → Fracture energy : 1.1982e+01 J
  → Total energy    : 2.3126e+02 J
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
  T_new: min=263.15 K, max=1023.18 K, mean=680.22 K
  T^n (self.T): min=263.15 K, max=1023.18 K
  ||ΔT||/||T|| = 1.447e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.310e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.307e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.600e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=680.22 K
  T^n (self.T): min=263.15 K, max=1023.18 K
  ||ΔT||/||T|| = 7.235e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.336e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.316e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.727e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=680.22 K
  T^n (self.T): min=263.15 K, max=1023.18 K
  ||ΔT||/||T|| = 3.617e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.015e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.773e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.160e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=680.22 K
  T^n (self.T): min=263.15 K, max=1023.18 K
  ||ΔT||/||T|| = 1.809e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.822e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.216e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.521e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1919e+02 J
  → Fracture energy : 1.2007e+01 J
  → Total energy    : 2.3120e+02 J
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
  T_new: min=263.15 K, max=1023.18 K, mean=679.36 K
  T^n (self.T): min=263.15 K, max=1023.18 K
  ||ΔT||/||T|| = 1.433e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.308e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.280e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.600e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=679.36 K
  T^n (self.T): min=263.15 K, max=1023.18 K
  ||ΔT||/||T|| = 7.167e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.332e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.291e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.727e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=679.36 K
  T^n (self.T): min=263.15 K, max=1023.18 K
  ||ΔT||/||T|| = 3.583e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.012e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.755e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.161e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=679.36 K
  T^n (self.T): min=263.15 K, max=1023.18 K
  ||ΔT||/||T|| = 1.792e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.794e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.203e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.523e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1911e+02 J
  → Fracture energy : 1.2031e+01 J
  → Total energy    : 2.3115e+02 J
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
  T_new: min=263.15 K, max=1023.18 K, mean=678.52 K
  T^n (self.T): min=263.15 K, max=1023.18 K
  ||ΔT||/||T|| = 1.420e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.305e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.250e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.581e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=678.52 K
  T^n (self.T): min=263.15 K, max=1023.18 K
  ||ΔT||/||T|| = 7.100e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.327e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.262e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.707e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=678.52 K
  T^n (self.T): min=263.15 K, max=1023.18 K
  ||ΔT||/||T|| = 3.550e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.008e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.734e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.147e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=678.52 K
  T^n (self.T): min=263.15 K, max=1023.18 K
  ||ΔT||/||T|| = 1.775e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.768e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.189e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.513e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1903e+02 J
  → Fracture energy : 1.2056e+01 J
  → Total energy    : 2.3109e+02 J
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
  T_new: min=263.15 K, max=1023.18 K, mean=677.68 K
  T^n (self.T): min=263.15 K, max=1023.18 K
  ||ΔT||/||T|| = 1.407e-03
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.302e-04
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.216e-03
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.545e-02

Convergence check


#### Iteration 2/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=677.68 K
  T^n (self.T): min=263.15 K, max=1023.18 K
  ||ΔT||/||T|| = 7.035e-05
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.323e-05
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 2.235e-04
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.670e-03

Convergence check


#### Iteration 3/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=677.68 K
  T^n (self.T): min=263.15 K, max=1023.18 K
  ||ΔT||/||T|| = 3.517e-06
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 1.004e-06
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.715e-05
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 2.118e-04

Convergence check


#### Iteration 4/200


**[INFO]** Assembling thermal problem...

  Building weak form, volume integrals (dx) for uo2, tag = 10
  → q_third[uo2](W/m3) min = 0.00e+00, max = 0.00e+00, mean = 0.00e+00
  Linear solver
  T_new: min=263.15 K, max=1023.18 K, mean=677.68 K
  T^n (self.T): min=263.15 K, max=1023.18 K
  ||ΔT||/||T|| = 1.759e-07
  [adaptive] relax_T=0.95

**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 40 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 30 → 0.0
  Building weak form, volume integrals (dx) for uo2, tag = 10
  Linear solver
  ||Δu||/||u|| = 6.742e-08
  [adaptive] relax_u=0.95

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'uo2' material
  - Material 'uo2': AT1 solve. Gc=3.72e+02, sigma_c=1.00e+09
  ||ΔD||/||D|| = 1.177e-06
  [adaptive] relax_D=0.95
  |ΔD|_∞ = 1.494e-05

Convergence check

**[SUCCESS]** Staggered solver converged in 4 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.1895e+02 J
  → Fracture energy : 1.2080e+01 J
  → Total energy    : 2.3103e+02 J
Exporting results to VTU file...
  → Projecting result fields for all materials...
  → Adding damage field to VTU...
VTU file exported to: output/fields_0099.vtu

Simulation completed in 508.07 s
Total time steps solved: 100
