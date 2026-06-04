Info    : Reading 'mesh.msh'...
Info    : 34 entities
Info    : 4645 nodes
Info    : 26744 elements
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
  → Time steps          : 50
  → Regime              : 3d
  → Models active       :
      thermal    → OFF
      mechanical → ON
      damage     → ON
      cluster    → OFF
      plasticity → OFF
  → Gap conductance     : None (value = 0.0)



### FiniteElementSetup initializer

Mechanical element order: 1
Thermal function space (V_t): FunctionSpace(<Mesh #0>, Basix element (P, tetrahedron, 1, gll_warped, unset, False, float64, []))
Mechanical function space (V_m): FunctionSpace(<Mesh #0>, blocked element (Basix element (P, tetrahedron, 1, gll_warped, unset, False, float64, []), (3,)))
Scalar function space (V_d): FunctionSpace(<Mesh #0>, Basix element (P, tetrahedron, 1, gll_warped, unset, False, float64, []))
Scalar function space (Q): FunctionSpace(<Mesh #0>, Basix element (P, tetrahedron, 0, gll_warped, unset, True, float64, []))
[Solver] initializer
  Applied relaxation factor:
  → Temperature  : 0.9
  → Displacement : 0.4
  → Damage       : 0.4
  Adaptive relaxation enabled
  → relax_growth  : 1.1
  → relax_shrink : 0.7
  → relax_min  : 0.05
  → relax_max : 1.0


[MechanicalModel] initializer
[MechanicalModel] options loaded from input.yaml:
  solver              : linear
  linear_solver       : iterative_amg
  rtol                : 1e-06
  stag_tol            : 1e-06
  convergence         : rel_norm
  debug               : False
DamageModel initializer
Options loaded from input.yaml:
  type                : AT1
  solver              : linear
  linear_solver       : iterative_amg
  lc                  : 5e-07
  hybrid_constraint   : True
  rtol                : 1e-06
  stag_tol            : 1e-06
  convergence         : rel_norm
[spine.load_materials]
Material loaded: solid
  → k defined as constant: 5.0
  → Gc defined as constant: 2.0
  - Material 'solid': sigma_c (AT1) from Gc = 2.00 J/m2
  → constitutive model: lame
  E               → 385000000000.0 (float)
  G               → 156504065040.65042 (float)
  Gc              → 2.0 (float)
  T_initial       → 1023.15 (float)
  T_ref           → 298.15 (float)
  alpha           → 1e-05 (float)
  bulk_modulus    → 237654320987.6543 (float)
  constitutive_mode → lame (str)
  cp              → 280.0 (float)
  k               → 5.0 (float)
  lmbda           → 133318277627.22072 (float)
  name            → UO2 (str)
  nu              → 0.23 (float)
  rho             → 10970.0 (float)
  sigma_c         → 759934207.6785332 (float)
[spine.initialize_fields]

Initializing the displacement field...
  Initial u: min=0.00e+00 m, max=0.00e+00 m, mean=0.00e+00 m

Initializing the damage field...



***


### spine - set_boundary_conditions


***



Loading boundary conditions from 'boundary_conditions.yaml'
  **[INFO]** Clamp_x mechanical BC on 'solid' → 0.0 (first step) at region 'xmin'
  **[INFO]** Clamp_y mechanical BC on 'solid' → 0.0 (first step) at region 'ymin'
  **[INFO]** Clamp_z mechanical BC on 'solid' → 0.0 (first step) at region 'zmin'
  **[INFO]** Dirichlet_z mechanical BC on 'solid' → 0.0 (first step) at region 'zmax'

Setting damage boundary conditions...
Computing symbolic result fields (strain, stress, ...)

**[INFO]** Hot-reload of allow-listed input.yaml parameters is active. Edit input.yaml during the run; changes apply at the next step boundary. Allowed keys: damage.{stag_tol,rtol,hybrid_constraint,gamma_star}, mechanical.{stag_tol,rtol}, thermal.{stag_tol,rtol}, solver_settings.{max_iters,relax_*}.


## Step 01/50: t = 0.00e+00 s | LHR = 0.00e+00 W/m

  → dt=0: solving static step / initial condition
Computing symbolic result fields (strain, stress, ...)



***


### spine - solve


***



Current step = 0 | dt = 0.00e+00 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 0.0
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=0.40

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.40
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 1 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 0.0000e+00 J
  → Fracture energy : 0.0000e+00 J
  → Total energy    : 0.0000e+00 J


## Step 02/50: t = 2.04e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 1 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.040816326530612e-08
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.000e+00
  [adaptive] relax_u=0.40

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.000e+00
  [adaptive] relax_D=0.40
  |ΔD|_∞ = 4.000e-01

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.040816326530612e-08
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.344e-01
  [adaptive] relax_u=0.44

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 5.657e-01
  [adaptive] relax_D=0.44
  |ΔD|_∞ = 3.155e-01

Convergence check


#### Iteration 3/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.040816326530612e-08
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.867e-01
  [adaptive] relax_u=0.48

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 6.475e-02
  [adaptive] relax_D=0.48
  |ΔD|_∞ = 1.584e-01

Convergence check


#### Iteration 4/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.040816326530612e-08
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.766e-01
  [adaptive] relax_u=0.53

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.53
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 5/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.040816326530612e-08
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.003e-01
  [adaptive] relax_u=0.59

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.59
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 6/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.040816326530612e-08
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.157e-02
  [adaptive] relax_u=0.64

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.64
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 7/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.040816326530612e-08
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.350e-02
  [adaptive] relax_u=0.71

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.71
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 8/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.040816326530612e-08
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.199e-03
  [adaptive] relax_u=0.78

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.78
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 9/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.040816326530612e-08
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.948e-03
  [adaptive] relax_u=0.86

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.86
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 10/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.040816326530612e-08
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.152e-04
  [adaptive] relax_u=0.94

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.94
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 11/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.040816326530612e-08
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.122e-04
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 12/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.040816326530612e-08
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.756e-06
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 13/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.040816326530612e-08
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.960e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 13 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 9.2132e-09 J
  → Fracture energy : 6.4857e-08 J
  → Total energy    : 7.4070e-08 J


## Step 03/50: t = 4.08e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 2 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.081632653061224e-08
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.737e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.081632653061224e-08
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.172e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.70
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.5949e-08 J
  → Fracture energy : 6.4857e-08 J
  → Total energy    : 1.0081e-07 J


## Step 04/50: t = 6.12e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 3 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.122448979591837e-08
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.333e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.70
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.122448979591837e-08
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.249e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.49
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.0886e-08 J
  → Fracture energy : 6.4857e-08 J
  → Total energy    : 1.4574e-07 J


## Step 05/50: t = 8.16e-02 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 4 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.163265306122448e-08
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.500e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.49
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.163265306122448e-08
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.183e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.34
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.4380e-07 J
  → Fracture energy : 6.4857e-08 J
  → Total energy    : 2.0865e-07 J


## Step 06/50: t = 1.02e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 5 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.020408163265306e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.000e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.34
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.020408163265306e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.414e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.24
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2468e-07 J
  → Fracture energy : 6.4857e-08 J
  → Total energy    : 2.8954e-07 J


## Step 07/50: t = 1.22e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 6 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2244897959183673e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.667e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.24
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.2244897959183673e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.267e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.17
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.2354e-07 J
  → Fracture energy : 6.4857e-08 J
  → Total energy    : 3.8840e-07 J


## Step 08/50: t = 1.43e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 7 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.4285714285714285e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.429e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.17
  |ΔD|_∞ = 0.000e+00

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.4285714285714285e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.407e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 0.000e+00
  [adaptive] relax_D=0.12
  |ΔD|_∞ = 0.000e+00

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.4038e-07 J
  → Fracture energy : 6.4857e-08 J
  → Total energy    : 5.0524e-07 J


## Step 09/50: t = 1.63e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 8 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.6326530612244896e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.250e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.170e-02
  [adaptive] relax_D=0.12
  |ΔD|_∞ = 2.669e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.6326530612244896e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.420e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.033e-02
  [adaptive] relax_D=0.13
  |ΔD|_∞ = 2.355e-02

Convergence check


#### Iteration 3/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.6326530612244896e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.213e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.002e-02
  [adaptive] relax_D=0.14
  |ΔD|_∞ = 2.286e-02

Convergence check


#### Iteration 4/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.6326530612244896e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.211e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 9.598e-03
  [adaptive] relax_D=0.16
  |ΔD|_∞ = 2.189e-02

Convergence check


#### Iteration 5/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.6326530612244896e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.791e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 9.055e-03
  [adaptive] relax_D=0.17
  |ΔD|_∞ = 2.065e-02

Convergence check


#### Iteration 6/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.6326530612244896e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.210e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 8.401e-03
  [adaptive] relax_D=0.19
  |ΔD|_∞ = 1.916e-02

Convergence check


#### Iteration 7/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.6326530612244896e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.213e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 7.649e-03
  [adaptive] relax_D=0.21
  |ΔD|_∞ = 1.745e-02

Convergence check


#### Iteration 8/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.6326530612244896e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.797e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 6.820e-03
  [adaptive] relax_D=0.23
  |ΔD|_∞ = 1.555e-02

Convergence check


#### Iteration 9/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.6326530612244896e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.210e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 5.938e-03
  [adaptive] relax_D=0.25
  |ΔD|_∞ = 1.354e-02

Convergence check


#### Iteration 10/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.6326530612244896e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.980e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 5.034e-03
  [adaptive] relax_D=0.28
  |ΔD|_∞ = 1.148e-02

Convergence check


#### Iteration 11/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.6326530612244896e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.007e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 4.141e-03
  [adaptive] relax_D=0.31
  |ΔD|_∞ = 9.445e-03

Convergence check


#### Iteration 12/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.6326530612244896e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.189e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 3.292e-03
  [adaptive] relax_D=0.34
  |ΔD|_∞ = 7.507e-03

Convergence check


#### Iteration 13/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.6326530612244896e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.181e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.516e-03
  [adaptive] relax_D=0.37
  |ΔD|_∞ = 5.738e-03

Convergence check


#### Iteration 14/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.6326530612244896e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.418e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.839e-03
  [adaptive] relax_D=0.41
  |ΔD|_∞ = 4.193e-03

Convergence check


#### Iteration 15/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.6326530612244896e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.184e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.276e-03
  [adaptive] relax_D=0.45
  |ΔD|_∞ = 2.909e-03

Convergence check


#### Iteration 16/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.6326530612244896e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.210e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 8.333e-04
  [adaptive] relax_D=0.49
  |ΔD|_∞ = 1.901e-03

Convergence check


#### Iteration 17/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.6326530612244896e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.836e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 5.071e-04
  [adaptive] relax_D=0.54
  |ΔD|_∞ = 1.157e-03

Convergence check


#### Iteration 18/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.6326530612244896e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.809e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.837e-04
  [adaptive] relax_D=0.59
  |ΔD|_∞ = 6.470e-04

Convergence check


#### Iteration 19/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.6326530612244896e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.058e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.434e-04
  [adaptive] relax_D=0.65
  |ΔD|_∞ = 3.270e-04

Convergence check


#### Iteration 20/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.6326530612244896e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 0.000e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 6.392e-05
  [adaptive] relax_D=0.72
  |ΔD|_∞ = 1.458e-04

Convergence check


#### Iteration 21/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.6326530612244896e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.212e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.432e-05
  [adaptive] relax_D=0.79
  |ΔD|_∞ = 5.547e-05

Convergence check


#### Iteration 22/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.6326530612244896e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.280e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 7.503e-06
  [adaptive] relax_D=0.87
  |ΔD|_∞ = 1.711e-05

Convergence check


#### Iteration 23/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.6326530612244896e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.349e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.721e-06
  [adaptive] relax_D=0.96
  |ΔD|_∞ = 3.925e-06

Convergence check


#### Iteration 24/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.6326530612244896e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.186e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.449e-07
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.586e-07

Convergence check

**[SUCCESS]** Staggered solver converged in 24 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 5.7411e-07 J
  → Fracture energy : 6.5951e-08 J
  → Total energy    : 6.4007e-07 J


## Step 10/50: t = 1.84e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 9 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.8367346938775507e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.096e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 4.183e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.812e-01

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1.8367346938775507e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 7.056e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.534e-15
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.109e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.1453e-07 J
  → Fracture energy : 7.7465e-08 J
  → Total energy    : 7.9199e-07 J


## Step 11/50: t = 2.04e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 10 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.040816326530612e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.673e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 4.744e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.507e-01

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.040816326530612e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.752e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.935e-15
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.272e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 8.1915e-07 J
  → Fracture energy : 1.2861e-07 J
  → Total energy    : 9.4776e-07 J


## Step 12/50: t = 2.24e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 11 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.2448979591836733e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.412e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 3.796e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.276e-01

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.2448979591836733e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.554e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.111e-15
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.995e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 7.6578e-07 J
  → Fracture energy : 2.4590e-07 J
  → Total energy    : 1.0117e-06 J


## Step 13/50: t = 2.45e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 12 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.4489795918367347e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.409e-01
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.172e-01
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.988e-01

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.4489795918367347e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.111e-16
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 5.562e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.273e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7004e-07 J
  → Fracture energy : 3.3817e-07 J
  → Total energy    : 6.0821e-07 J


## Step 14/50: t = 2.65e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 13 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.6530612244897955e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.048e+00
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.420e-02
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.926e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.6530612244897955e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.415e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.257e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.7152e-08 J
  → Fracture energy : 3.4222e-07 J
  → Total energy    : 3.5937e-07 J


## Step 15/50: t = 2.86e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 14 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.857142857142857e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.547e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 5.102e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.146e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 2.857142857142857e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.022e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.036e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 1.9601e-08 J
  → Fracture energy : 3.4224e-07 J
  → Total energy    : 3.6184e-07 J


## Step 16/50: t = 3.06e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 15 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.061224489795918e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.684e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 3.154e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.470e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.061224489795918e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.601e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 9.592e-17
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.882e-16

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.2422e-08 J
  → Fracture energy : 3.4227e-07 J
  → Total energy    : 3.6469e-07 J


## Step 17/50: t = 3.27e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 16 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.265306122448979e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.258e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 9.257e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.025e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.265306122448979e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.813e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.878e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.998e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.5092e-08 J
  → Fracture energy : 3.4239e-07 J
  → Total energy    : 3.6749e-07 J


## Step 18/50: t = 3.47e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 17 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.4693877551020406e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.915e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.227e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.408e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.4693877551020406e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.548e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.190e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 2.7590e-08 J
  → Fracture energy : 3.4258e-07 J
  → Total energy    : 3.7017e-07 J


## Step 19/50: t = 3.67e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 18 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.6734693877551015e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.631e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.167e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.035e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.6734693877551015e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.086e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.150e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.0002e-08 J
  → Fracture energy : 3.4279e-07 J
  → Total energy    : 3.7279e-07 J


## Step 20/50: t = 3.88e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 19 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.877551020408163e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.337e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.430e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.165e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 3.877551020408163e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.864e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.081e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.332e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.1931e-08 J
  → Fracture energy : 3.4309e-07 J
  → Total energy    : 3.7502e-07 J


## Step 21/50: t = 4.08e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 20 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.081632653061224e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.071e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.093e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.000e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.081632653061224e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.051e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.470e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.998e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.4157e-08 J
  → Fracture energy : 3.4333e-07 J
  → Total energy    : 3.7748e-07 J


## Step 22/50: t = 4.29e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 21 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.285714285714285e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.819e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.709e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.347e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.285714285714285e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.089e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.341e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.887e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.5766e-08 J
  → Fracture energy : 3.4367e-07 J
  → Total energy    : 3.7944e-07 J


## Step 23/50: t = 4.49e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 22 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.4897959183673465e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.615e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.608e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.430e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.4897959183673465e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.755e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.494e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.7177e-08 J
  → Fracture energy : 3.4402e-07 J
  → Total energy    : 3.8119e-07 J


## Step 24/50: t = 4.69e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 23 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.6938775510204074e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.405e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.572e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.839e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.6938775510204074e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.090e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 3.022e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.8575e-08 J
  → Fracture energy : 3.4435e-07 J
  → Total energy    : 3.8292e-07 J


## Step 25/50: t = 4.90e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 24 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.897959183673469e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.212e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.645e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.811e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 4.897959183673469e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.087e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.755e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 3.9699e-08 J
  → Fracture energy : 3.4470e-07 J
  → Total energy    : 3.8440e-07 J


## Step 26/50: t = 5.10e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 25 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.10204081632653e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.042e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.787e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.582e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.10204081632653e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.315e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.875e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.776e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.0348e-08 J
  → Fracture energy : 3.4509e-07 J
  → Total energy    : 3.8544e-07 J


## Step 27/50: t = 5.31e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 26 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.306122448979591e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.884e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 9.787e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.107e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.306122448979591e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 4.295e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.794e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.1694e-08 J
  → Fracture energy : 3.4535e-07 J
  → Total energy    : 3.8704e-07 J


## Step 28/50: t = 5.51e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 27 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.510204081632652e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.731e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.319e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.964e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.510204081632652e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 9.592e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.505e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.0704e-08 J
  → Fracture energy : 3.4588e-07 J
  → Total energy    : 3.8658e-07 J


## Step 29/50: t = 5.71e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 28 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.714285714285714e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.607e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.799e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.047e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.714285714285714e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.655e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.269e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.3223e-08 J
  → Fracture energy : 3.4595e-07 J
  → Total energy    : 3.8917e-07 J


## Step 30/50: t = 5.92e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 29 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.918367346938775e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.457e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.560e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.170e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 5.918367346938775e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.742e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.631e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.0493e-08 J
  → Fracture energy : 3.4662e-07 J
  → Total energy    : 3.8711e-07 J


## Step 31/50: t = 6.12e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 30 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.122448979591836e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.353e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 4.309e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 5.059e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.122448979591836e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.949e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.492e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.443e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.3298e-08 J
  → Fracture energy : 3.4662e-07 J
  → Total energy    : 3.8992e-07 J


## Step 32/50: t = 6.33e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 31 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.326530612244898e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.229e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.474e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.115e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.326530612244898e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.034e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 3.134e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.220e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.2069e-08 J
  → Fracture energy : 3.4707e-07 J
  → Total energy    : 3.8914e-07 J


## Step 33/50: t = 6.53e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 32 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.530612244897958e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.136e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.610e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.640e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.530612244897958e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.329e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.857e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.331e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.4344e-08 J
  → Fracture energy : 3.4712e-07 J
  → Total energy    : 3.9147e-07 J


## Step 34/50: t = 6.73e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 33 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.734693877551019e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.033e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.405e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.907e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.734693877551019e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.369e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.997e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.1056e-08 J
  → Fracture energy : 3.4773e-07 J
  → Total energy    : 3.8879e-07 J


## Step 35/50: t = 6.94e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 34 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.938775510204081e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.948e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 3.354e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.451e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 6.938775510204081e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.224e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.636e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.998e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.3560e-08 J
  → Fracture energy : 3.4774e-07 J
  → Total energy    : 3.9130e-07 J


## Step 36/50: t = 7.14e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 35 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.142857142857142e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.860e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 7.057e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 3.751e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.142857142857142e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.569e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 3.421e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.887e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.3818e-08 J
  → Fracture energy : 3.4796e-07 J
  → Total energy    : 3.9178e-07 J


## Step 37/50: t = 7.35e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 36 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.346938775510203e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.785e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.267e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.175e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.346938775510203e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 8.935e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 3.333e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.331e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.3141e-08 J
  → Fracture energy : 3.4826e-07 J
  → Total energy    : 3.9140e-07 J


## Step 38/50: t = 7.55e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 37 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.551020408163265e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.705e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.883e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.856e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.551020408163265e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.567e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.593e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.109e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.4689e-08 J
  → Fracture energy : 3.4834e-07 J
  → Total energy    : 3.9303e-07 J


## Step 39/50: t = 7.76e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 38 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.755102040816326e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.634e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.583e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.282e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.755102040816326e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.177e-18
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.744e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.2582e-08 J
  → Fracture energy : 3.4874e-07 J
  → Total energy    : 3.9132e-07 J


## Step 40/50: t = 7.96e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 39 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.959183673469387e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.565e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 4.374e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.129e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 7.959183673469387e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.256e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.161e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.776e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.4769e-08 J
  → Fracture energy : 3.4874e-07 J
  → Total energy    : 3.9351e-07 J


## Step 41/50: t = 8.16e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 40 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.163265306122447e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.504e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.128e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 7.061e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.163265306122447e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.436e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 3.354e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.887e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.3153e-08 J
  → Fracture energy : 3.4907e-07 J
  → Total energy    : 3.9222e-07 J


## Step 42/50: t = 8.37e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 41 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.367346938775509e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.440e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 4.227e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.271e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.367346938775509e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.935e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.528e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.109e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.5217e-08 J
  → Fracture energy : 3.4908e-07 J
  → Total energy    : 3.9430e-07 J


## Step 43/50: t = 8.57e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 42 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.57142857142857e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.383e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.650e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.268e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.57142857142857e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 6.553e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.867e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.887e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.2257e-08 J
  → Fracture energy : 3.4948e-07 J
  → Total energy    : 3.9174e-07 J


## Step 44/50: t = 8.78e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 43 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.775510204081631e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.325e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.746e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.504e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.775510204081631e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.231e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.339e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.4281e-08 J
  → Fracture energy : 3.4948e-07 J
  → Total energy    : 3.9377e-07 J


## Step 45/50: t = 8.98e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 44 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.979591836734693e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.275e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 8.722e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.148e-02

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 8.979591836734693e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.108e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 3.249e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.998e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.3505e-08 J
  → Fracture energy : 3.4969e-07 J
  → Total energy    : 3.9319e-07 J


## Step 46/50: t = 9.18e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 45 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.183673469387754e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.226e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 4.855e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 4.327e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.183673469387754e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.192e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 3.074e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.109e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.4065e-08 J
  → Fracture energy : 3.4979e-07 J
  → Total energy    : 3.9386e-07 J


## Step 47/50: t = 9.39e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 46 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.387755102040815e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.174e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 8.389e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 6.761e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.387755102040815e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.386e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 3.899e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.554e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.3105e-08 J
  → Fracture energy : 3.5000e-07 J
  → Total energy    : 3.9310e-07 J


## Step 48/50: t = 9.59e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 47 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.591836734693876e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.127e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.223e-04
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.687e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.591836734693876e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 1.272e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 3.173e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.665e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.4181e-08 J
  → Fracture energy : 3.5006e-07 J
  → Total energy    : 3.9424e-07 J


## Step 49/50: t = 9.80e-01 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 48 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.795918367346939e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.084e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 1.117e-03
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 8.534e-03

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 9.795918367346939e-07
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 3.803e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 3.379e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 2.442e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.2124e-08 J
  → Fracture energy : 3.5032e-07 J
  → Total energy    : 3.9244e-07 J


## Step 50/50: t = 1.00e+00 s | LHR = 0.00e+00 W/m




***


### spine - solve


***



Current step = 49 | dt = 2.04e-02 s
Coupling = staggered
  → Max iterations              : 100
  → Staggering tolerance |ΔT|   : 1.0e-04
  → Staggering tolerance |Δu|   : 1.0e-06
  → Staggering tolerance |ΔD|   : 1.0e-06
  → Relative tolerance th       : 1.0e-06
  → Relative tolerance mech     : 1.0e-06
  → Relative tolerance dmg      : 1.0e-06


#### Iteration 1/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1e-06
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 2.040e-02
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.302e-05
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.999e-04

Convergence check


#### Iteration 2/100


**[INFO]** Assembling mechanical problem...
  **[INFO]** Updating Displacement Dirichlet on region 3 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 4 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 7 → 0.0
  **[INFO]** Updating Displacement Dirichlet on region 5 → 1e-06
  Building weak form, volume integrals (dx) for solid, tag = 1
  Linear solver
  ||Δu||/||u|| = 5.031e-17
  [adaptive] relax_u=1.00

**[INFO]** Assembling damage (AT1) problem...
Solving damage problem for 'solid' material
  - Material 'solid': AT1 solve. Gc=2.00e+00, sigma_c=7.60e+08
  ||ΔD||/||D|| = 2.541e-16
  [adaptive] relax_D=1.00
  |ΔD|_∞ = 1.887e-15

Convergence check

**[SUCCESS]** Staggered solver converged in 2 iterations.
Computing symbolic result fields (strain, stress, ...)
  → Elastic energy  : 4.3851e-08 J
  → Fracture energy : 3.5032e-07 J
  → Total energy    : 3.9417e-07 J

Simulation completed in 92.85 s
Total time steps solved: 50
